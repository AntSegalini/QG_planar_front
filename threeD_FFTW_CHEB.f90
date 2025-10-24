! thanks to https://stackoverflow.com/questions/41730864/incorrect-output-via-fftw-mpi-r2c-2d-and-fftw-mpi-c2r-2d

module threeD_FFTW_CHEB    
    ! This module creates the routines to perform the spectral transform of a 3D field using FFTW and Chebyshev polynomials
    ! The grid field has dimensions (Nz,Nx,Ny) and the spectral field has dimensions (Nz,Ny,Nx/2) (no alias correction) or (2Nz/3,2Ny/3, Nx/3) (with alias correction)
    ! The MPI parallelisation is done in the y direction (grid space) or in kx direction (spectral space)

    use mpi
    use, intrinsic :: iso_c_binding
    implicit none    
    include 'fftw3-mpi.f03'
    !---------------------------------------------------------------
    
    public ::   initialize_transform, terminate_FFTW, grid_2_spectral, get_x, get_y, get_z, get_spectral_shape, &
                get_kx, get_ky, grid_2_spectral_FFT, spectral_2_grid_FFT, spectral_2_grid, get_matrices, &
                CHEB_first_derivative_spectral, CHEB_second_derivative_spectral, get_grid_shape, grid_2_spectral_XYplane, &
                spectral_2_grid_XYplane

    private

    logical, save :: dealias
    integer, save :: local_Nx, local_x_offset, local_Ny, local_y_offset, Nx, Ny, Nz, Nx_t, Ny_t, Nz_t, kxm, kym, kzm
    integer, dimension(:), allocatable, save :: kx, ky, kx_trunc_I_OK, ky_trunc_I_OK
    real(8), parameter :: pi=4*atan(1.0_8)
    real(8), dimension(:,:), allocatable, save :: T, Tinv

    ! FFTW variables
    type(c_ptr), save :: plan_forward, plan_backward, plan_CHEBY
    real(c_double), dimension(:,:), pointer, save :: data
    complex(c_double_complex), dimension(:,:), pointer, save :: data_hat
    real(c_double), dimension(:,:,:), allocatable, save :: data2

    ! MPI variables
    integer, save :: ierr , rank , nproc

    contains

    ! #######################################################

    subroutine initialize_transform(Nx_input, Ny_input, Nz_input, dealiasing, ICOM_input)
        
        logical, intent(in) :: dealiasing
        integer, intent(in) :: Nx_input, Ny_input, Nz_input
        integer :: i, j, ICOM_input, info
        integer, dimension(:), allocatable :: ipiv      ! Pivot indices
        real(8), dimension(:), allocatable :: work, z   ! Workspace
        
        integer(C_INTPTR_T) :: alloc_local, local_Nx_c, local_x_offset_c, local_Ny_c, local_y_offset_c
        type(c_ptr) :: cdata_hat, cdata 

        ! Initialize MPI
        call MPI_Comm_rank(ICOM_input, rank, ierr)
        call MPI_Comm_size(ICOM_input, nproc, ierr)
        call fftw_mpi_init()

        ! ########################################################

        Nx=Nx_input;        Ny=Ny_input;        Nz=Nz_input;           dealias=dealiasing

        ! ########################################################

        ! Fourier part (2D MPI)
        alloc_local = fftw_mpi_local_size_2d(int(Nx/2+1,kind=C_INTPTR_T), int(Ny,kind=C_INTPTR_T), ICOM_input, &
                                            local_Nx_c, local_x_offset_c)
        cdata_hat = fftw_alloc_complex(alloc_local)
        call c_f_pointer(cdata_hat, data_hat, [int(Ny,kind=C_INTPTR_T),local_Nx_c])

        alloc_local = fftw_mpi_local_size_2d(int(Ny,kind=C_INTPTR_T), int(Nx/2+1,kind=C_INTPTR_T), ICOM_input, &
                                            local_Ny_c, local_y_offset_c)
        cdata = fftw_alloc_real(2*alloc_local)
        call c_f_pointer(cdata, data, [int(2*(Nx/2+1),kind=C_INTPTR_T),local_Ny_c])

        plan_forward = fftw_mpi_plan_dft_r2c_2d(int(Ny,kind=C_INTPTR_T), int(Nx,kind=C_INTPTR_T), data, data_hat, &
                                                ICOM_input, ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT))
        plan_backward = fftw_mpi_plan_dft_c2r_2d(int(Ny,kind=C_INTPTR_T), int(Nx,kind=C_INTPTR_T), data_hat, data, &
                                                ICOM_input, ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN))
        
        local_Nx = local_Nx_c                ! kx direction
        local_x_offset = local_x_offset_c
        local_Ny = local_Ny_c                ! y direction
        local_y_offset = local_y_offset_c

        ! wavenumbers and truncation
        kx=[(local_x_offset+j,j=0,local_Nx-1)]
        ky=[(j,j=0,Ny-1)]    
        where (ky>=Ny/2)        ky=ky-Ny
        if (dealias) then
            kx_trunc_I_OK   = pack([(j,j=1,local_Nx)],kx<Nx/3)        ! indices where kx is OK
            ky_trunc_I_OK   = pack([(j,j=1,Ny)],abs(ky)<Ny/3)         ! indices where ky is OK
            Nz_t=(2*Nz)/3
        else
            kx_trunc_I_OK   = pack([(j,j=1,local_Nx)],kx<Nx/2)        ! indices where kx is OK
            ky_trunc_I_OK   = pack([(j,j=1,Ny)],abs(ky)<Ny/2)         ! indices where ky is OK
            Nz_t=Nz
        end if
        Nx_t=size(kx_trunc_I_OK)
        Ny_t=size(ky_trunc_I_OK)

        ! ##############################################################################

        ! Gauss-Lobatto grid (Chebyshev)
        allocate(data2(Nz,Nx,local_Ny))
        plan_CHEBY = fftw_plan_many_r2r(1, (/Nz/), Nx*local_Ny, data2, (/Nz/), 1, Nz, &
                                        data2, (/Nz/), 1, Nz, (/FFTW_REDFT00/) , FFTW_MEASURE)


        ! Chebyshev polynomials (for the full discretization , T converts from spectral to grid)
        allocate(T(Nz,Nz),Tinv(Nz,Nz))
        z=get_z(Nz)
        T(:,1)=1
        T(:,2)=z
        do i=3,Nz
            T(:,i) = 2*z*T(:,i-1)  - T(:,i-2)
        end do

        Tinv=T ! inverse of T (Tinv converts from grid to spectral)
        allocate(ipiv(Nz),work(Nz))
        call dgetrf(Nz, Nz, Tinv, Nz, ipiv, info)
        call dgetri(Nz, Tinv, Nz, ipiv, work, Nz, info)
        deallocate(ipiv,work)

    end subroutine initialize_transform

    ! #######################################################

    subroutine get_matrices(T_t, T1_t, T2_t, T3_t, T4_t, Tinv_t)
        ! This subroutine computes the Chebyshev polynomials and their derivatives for the truncated discretization
        integer :: i, j, info
        integer, dimension(Nz_t) :: ipiv      ! Pivot indices
        real(8), dimension(Nz_t,Nz_t) :: T_t, T1_t, T2_t, T3_t, T4_t, Tinv_t    
        real(8), dimension(Nz_t) :: work, z      ! Workspace
        
        z=cos(pi*[(i,i=0,Nz_t-1)]/(Nz_t-1))

        ! Chebyshev polynomials
        T_t(:,1)=1
        T_t(:,2)=z
        T1_t(:,1)=0
        T1_t(:,2)=1
        T2_t(:,:2)=0
        T3_t(:,:2)=0
        T4_t(:,:2)=0
        do i=3,Nz_t
            T_t(:,i) = 2*z*T_t(:,i-1)  - T_t(:,i-2)
            T1_t(:,i)= 2*z*T1_t(:,i-1) + 2*T_t(:,i-1)  - T1_t(:,i-2)
            T2_t(:,i)= 2*z*T2_t(:,i-1) + 4*T1_t(:,i-1) - T2_t(:,i-2)
            T3_t(:,i)= 2*z*T3_t(:,i-1) + 6*T2_t(:,i-1) - T3_t(:,i-2)
            T4_t(:,i)= 2*z*T4_t(:,i-1) + 8*T3_t(:,i-1) - T4_t(:,i-2)
        end do

        ! inverse of T
        Tinv_t=T_t 
        call dgetrf(Nz_t, Nz_t, Tinv_t, Nz_t, ipiv, info)
        call dgetri(Nz_t, Tinv_t, Nz_t, ipiv, work, Nz_t, info)
        
    end subroutine get_matrices

    ! #######################################################

    function get_x()
        integer :: i
        real(8), dimension(Nx) :: get_x
        get_x = 2*pi*[(i,i=0,Nx-1)]/Nx
    end function get_x

    ! #######################################################

    function get_y()
        integer :: i
        real(8), dimension(local_Ny) :: get_y
        get_y = 2*pi*[(i+local_y_offset,i=0,local_Ny-1)]/Ny
    end function get_y
    
    ! #######################################################

    function get_z(NN)
        ! Chebyshev grid (full discretization)
        integer :: i,NN
        real(8), allocatable, dimension(:) :: get_z
        allocate(get_z(NN))
        get_z = cos(pi*[(i,i=0,NN-1)]/(NN-1))
    end function get_z
    
    ! #######################################################

    subroutine get_grid_shape(NGz,NGx,NGy)
        ! This subroutine returns the shape of the spectral field
        integer :: NGz,NGy,NGx
        NGz=Nz
        NGx=Nx
        NGy=local_Ny
    end subroutine

    ! #######################################################

    subroutine get_spectral_shape(Nzz,Nyy,Nxx)
        ! This subroutine returns the shape of the spectral field
        integer :: Nxx, Nyy, Nzz
        Nxx=Nx_t
        Nyy=Ny_t
        Nzz=Nz_t
    end subroutine

    ! #######################################################

    function get_kx()
        ! This function returns the wavenumbers in the x direction (MPI splitting)
        integer, dimension(Nx_t) :: get_kx
        get_kx = kx(kx_trunc_I_OK)
    end function get_kx

    ! #######################################################

    function get_ky()
        ! This function returns the wavenumbers in the y direction
        integer, dimension(Ny_t) :: get_ky
        get_ky = ky(ky_trunc_I_OK)
    end function get_ky

    ! #######################################################

    subroutine grid_2_spectral_XYplane(GRID, SPECTRAL)
        ! Routine that computes the transform of a grid to spectral space with the matrix method in z
        real(8), intent(in) :: GRID(:,:)
        complex(8), intent(inout) :: SPECTRAL(:,:)
        
        data(:Nx,:)=GRID      
        call fftw_mpi_execute_dft_r2c(plan_forward, data, data_hat)
        SPECTRAL=data_hat(ky_trunc_I_OK,:Nx_t)                  
    end subroutine 

    ! #######################################################
    
    subroutine spectral_2_grid_XYplane(SPECTRAL, GRID)
        ! Routine that computes the transform of spectral to grid space
        real(8), intent(inout) :: GRID(:,:)
        complex(8), intent(in) :: SPECTRAL(:,:)

        data_hat=0
        data_hat(ky_trunc_I_OK,:Nx_t)=SPECTRAL
        call fftw_mpi_execute_dft_c2r(plan_backward, data_hat, data)
        GRID=data(:Nx,:)/(Nx*Ny)        
    end subroutine 
    
    ! #######################################################
    
    subroutine grid_2_spectral(GRID, SPECTRAL)
        ! Routine that computes the transform of a grid to spectral space with the matrix method in z
        integer :: i, j
        real(8), intent(in) :: GRID(:,:,:)
        complex(8), intent(inout) :: SPECTRAL(:,:,:)
        
        do i=1,Nz_t
            ! Chebyshev transform (3x faster than matmul for each coordinate point)
            data(:Nx,:)=Tinv(i,1)*GRID(1,:,:)
            do j=2,Nz
                data(:Nx,:)=data(:Nx,:)+Tinv(i,j)*GRID(j,:,:)
            end do

            ! Fourier planar transform        
            call fftw_mpi_execute_dft_r2c(plan_forward, data, data_hat)
            SPECTRAL(i,:,:)=data_hat(ky_trunc_I_OK,:Nx_t)
        end do          
    end subroutine 

    ! #######################################################

    subroutine grid_2_spectral_FFT(GRID, SPECTRAL)
        ! Routine that computes the transform of a grid to spectral space with the FFT method in z
        integer :: i, j
        real(8), intent(in) :: GRID(:,:,:)
        complex(8), intent(inout) :: SPECTRAL(:,:,:)
        ! #######################################################
        ! Chebyshev transform
        data2=GRID
        call fftw_execute_r2r(plan_CHEBY,data2,data2)

        ! Fourier planar transform
        do i =1,Nz_t
            data(:Nx,:)=data2(i,:,:)
            call fftw_mpi_execute_dft_r2c(plan_forward, data, data_hat)
            SPECTRAL(i,:,:)=data_hat(ky_trunc_I_OK,:Nx_t)
        end do  
        
        SPECTRAL(1,:,:)=SPECTRAL(1,:,:)/2
        if (.not.dealias) SPECTRAL(Nz,:,:)=SPECTRAL(Nz,:,:)/2
        SPECTRAL=SPECTRAL/(Nz-1)        
    end subroutine 

    ! #######################################################

    subroutine spectral_2_grid(SPECTRAL, GRID)
        ! Routine that computes the transform of spectral to grid space with the matrix method in z
        integer :: i, j
        real(8), intent(inout) :: GRID(:,:,:)
        complex(8), intent(in) :: SPECTRAL(:,:,:)
        ! ####################################################### 

        ! Fourier planar transform
        data_hat=0
        do i =1,Nz_t
            data_hat(ky_trunc_I_OK,:Nx_t)=SPECTRAL(i,:,:)
            call fftw_mpi_execute_dft_c2r(plan_backward, data_hat, data)
            data2(i,:,:)=data(:Nx,:)/(Nx*Ny)  
        end do  
        
        do i=1,Nz
            GRID(i,:,:)=T(i,1)*data2(1,:,:)
            do j=2,Nz_t
                GRID(i,:,:)=GRID(i,:,:)+T(i,j)*data2(j,:,:)
            end do
        end do      
    end subroutine 
    
    ! #######################################################

    subroutine spectral_2_grid_FFT(SPECTRAL, GRID)
        ! Routine that computes the transform of spectral to grid space with the FFT method in z
        integer :: i, j
        real(8), intent(inout) :: GRID(:,:,:)
        complex(8), intent(in) :: SPECTRAL(:,:,:)
        ! ####################################################### 
        ! Fourier planar transform
        data_hat=0
        do i =1,Nz_t
            data_hat(ky_trunc_I_OK,:Nx_t)=SPECTRAL(i,:,:)
            call fftw_mpi_execute_dft_c2r(plan_backward, data_hat, data)
            data2(i,:,:)=data(:Nx,:)/(Nx*Ny)  
        end do  
        
        if (dealias) then
            data2(2:Nz_t,:,:)=data2(2:Nz_t,:,:)/2
            data2(Nz_t+1:,:,:)=0
        else
            data2(2:Nz-1,:,:)=data2(2:Nz-1,:,:)/2
        end if
        
        ! inverse Chebyshev transform
        call fftw_execute_r2r(plan_CHEBY,data2,data2) 
        GRID=data2   
    end subroutine 
        
    ! #######################################################

    subroutine CHEB_first_derivative_spectral(SPECTRAL,SPECTRAL_derivative)
        integer :: p, k
        complex(8) , intent(in) :: SPECTRAL(:,:,:)
        complex(8) , intent(inout) :: SPECTRAL_derivative(:,:,:)
        
        do k=1,Nz_t-1
            SPECTRAL_derivative(k+1,:,:)=k*SPECTRAL(k+1,:,:)
        end do 
        SPECTRAL_derivative(1,:,:)=sum(SPECTRAL_derivative(2::2,:,:),dim=1)
        do k=1,Nz_t-1
            SPECTRAL_derivative(k+1,:,:)=2*sum(SPECTRAL_derivative(k+2::2,:,:),dim=1)
        end do 
        
        if (.not.dealias) SPECTRAL_derivative(Nz,:,:) = 0
    end subroutine

    ! #######################################################

    subroutine CHEB_second_derivative_spectral(SPECTRAL,SPECTRAL_derivative)
        integer :: p, k
        complex(8) , intent(in) :: SPECTRAL(:,:,:)
        complex(8) , intent(inout) :: SPECTRAL_derivative(:,:,:)
        complex(8) , dimension(Nz_t,Ny_t,Nx_t) :: AUX
        
        do k=2,Nz_t-1
            SPECTRAL_derivative(k+1,:,:)=k**3*SPECTRAL(k+1,:,:)
            AUX(k+1,:,:)=k*SPECTRAL(k+1,:,:)
        end do 
        SPECTRAL_derivative(1,:,:)=sum(SPECTRAL_derivative(3::2,:,:),dim=1)/2
        do k=1,Nz_t-1
            SPECTRAL_derivative(k+1,:,:)=sum(SPECTRAL_derivative(k+3::2,:,:),dim=1)-k**2*sum(AUX(k+3::2,:,:),dim=1)
        end do 
        
        if (.not.dealias) SPECTRAL_derivative(Nz-1:,:,:) = 0
    end subroutine

    ! #######################################################

    subroutine terminate_FFTW()
        call fftw_destroy_plan(plan_forward)
        call fftw_destroy_plan(plan_backward)
        call fftw_destroy_plan(plan_CHEBY)
        call fftw_mpi_cleanup()
        ! deallocate(data, data_hat)
    end subroutine terminate_FFTW
 
end module threeD_FFTW_CHEB
