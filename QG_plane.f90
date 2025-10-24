program QG_plane

    use mpi
    use threeD_FFTW_CHEB
    use, intrinsic :: iso_c_binding        
    implicit none

    real(8), parameter :: pi = 4*atan(1.0_8)
    complex(8), parameter :: Iunit = (0.0_8, 1.0_8)

    ! ###################################################################

    logical :: dealiasing_condition, nonlinear_condition

    integer :: i, j, p, iter, mm, nn, &
               Nx, Ny, Nz, NSx, NSy, NSz, NGx, NGy, NGz, &
               ierr, rank, nproc, kx_initial, ky_initial, Nsave
    
    real(8), allocatable, dimension(:)   :: x, y, z, kx, ky
    real(8), allocatable, dimension(:,:) :: T, T1, T2, T3, T4, Tinv, phi, k2, OPE1, OPE2
    
    real(8) :: Lx, Ly, Lz, f0, beta, DT, nu, dPSI0bottom, dPSI0top, tempo, tempo1, tempo2, Tfinal, &
               C0, C, Az, Bz, Initial_perturbation, maxU, maxV, friction_bottom

    real(8), allocatable, dimension(:) :: detadz, Ubar, dUbar, d2Ubar, rhobar, drhobar, N2bar, dN2bar, Q0, & 
                                z_file, U0_file, dU0_file, d2U0_file, rho_file, drho_file, N2_file, dN2_file, &
                                Q0_extended, Ubar_extended, z_extended
                   
    complex(8), allocatable, dimension(:,:,:) :: Q_hat, Q_hat_old
    complex(8), allocatable, dimension(:,:)   :: dPSI_bottom_hat, dPSI_top_hat, dPSI_bottom_hat_old, &
                                                 dPSI_top_hat_old, A0, Ak2

    character(len=100) :: line, str2, fileName_root

    ! ###################################################################################

    ! MPI initialization
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! ###################################################################################

    ! reading the input parameters
    open (unit=10, file='INPUT.txt', status='old')
    read(10,*) ! skip the header
    read(10,*)  fileName_root, Nx, Ny, Nz, Lx, Ly, Lz, DT, f0, beta, nu, p, Initial_perturbation, &
                kx_initial, ky_initial, Nsave, Tfinal, dealiasing_condition, nonlinear_condition, friction_bottom
    close (10) 
    
    ! ###################################################################################

    ! create the output directory
    if (rank==0) then
        call system('mkdir -p '//adjustl(trim(fileName_root)))
        call system('cp INPUT.txt ./'//adjustl(trim(fileName_root))//'/')
        call system('cp PROFILES.txt ./'//adjustl(trim(fileName_root))//'/')
    end if

    ! ###################################################################################
    
    ! Initialize the transform
    call initialize_transform(Nx, Ny, Nz, dealiasing_condition, MPI_COMM_WORLD)
    call get_spectral_shape(NSz, NSy, NSx)
    call get_grid_shape(NGz, NGx, NGy)
    
    ! ###################################################################################

    ! Allocate the matrices for the operators and the coordinates
    allocate( T(NSz, NSz), T1(NSz, NSz), T2(NSz, NSz), T3(NSz, NSz), T4(NSz, NSz), Tinv(NSz, NSz))
    call get_matrices(T, T1, T2, T3, T4, Tinv)

    allocate( z_extended(NGz), Ubar_extended(NGz), Q0_extended(NGz) )
    x  = get_x() * Lx / (2*pi)
    y  = get_y() * Ly / (2*pi)
    kx = get_kx()* 2*pi / Lx
    ky = get_ky()* 2*pi / Ly
    z  = (1-get_z(NSz)) / 2 * Lz 
    z_extended = (1-get_z(NGz)) / 2 * Lz

    T1=T1*(-2/Lz) ! Transformation of the derivative operators
    T2=T2*(4/Lz**2)

    ! ###################################################################################

    ! useful operators
    allocate(k2(NSy,NSx), OPE1(NSy,NSx), OPE2(NSy,NSx))
    do i=1,NSy
        k2(i,:) = ky(i)**2
    end do
    do j=1,NSx
        k2(:,j) = k2(:,j) + kx(j)**2
    end do
    OPE1 = 1 - nu*DT*k2**p
    OPE2 = 1 + nu*DT*k2**p

    ! ###################################################################################
    
    ! read the U, rho, N2 profile file and interpolate
    nn = 0
    open (unit=10, file='PROFILES.txt', status='old', action='read')
    read(10, "(A)", iostat=mm) line ! Skip the header
    do
        read(10, "(A)", iostat=mm) line
        if (mm /= 0) exit  ! End of file or error
        nn = nn + 1
    end do

    allocate(z_file(nn), U0_file(nn), dU0_file(nn), d2U0_file(nn), rho_file(nn), drho_file(nn), N2_file(nn), dN2_file(nn))
    rewind(10)
    read(10, "(A)", iostat=mm) line ! Skip the header
    do i=1,nn
        read(10, *,iostat=mm) z_file(i), U0_file(i), dU0_file(i), d2U0_file(i), rho_file(i), drho_file(i), N2_file(i), dN2_file(i)
    end do
    close(10)

    allocate(Ubar(NSz), dUbar(NSz), d2Ubar(NSz), rhobar(NSz), drhobar(NSz), N2bar(NSz), dN2bar(NSz))
    do i=1,NSz
        Ubar(i)     = interp1d(z_file, U0_file,   z(i))
        dUbar(i)    = interp1d(z_file, dU0_file,  z(i))
        d2Ubar(i)   = interp1d(z_file, d2U0_file, z(i))
        rhobar(i)   = interp1d(z_file, rho_file,  z(i))
        drhobar(i)  = interp1d(z_file, drho_file, z(i))
        N2bar(i)    = interp1d(z_file, N2_file,   z(i))
        dN2bar(i)   = interp1d(z_file, dN2_file,  z(i))
    end do
    dN2_file = - f0**2/N2_file * ( (drho_file/rho_file-dN2_file/N2_file)*dU0_file + d2U0_file ) + beta
    do i=1,NGz
        Ubar_extended(i) = interp1d(z_file, U0_file,   z_extended(i))
        Q0_extended(i)   = interp1d(z_file, dN2_file,   z_extended(i))
    end do
    deallocate(z_file, U0_file, dU0_file, d2U0_file, rho_file, drho_file, N2_file, dN2_file)

    ! ###################################################################################

    ! background operators
    allocate(Q0(NSz))
    Q0 = - f0**2/N2bar *( (drhobar/rhobar-dN2bar/N2bar)*dUbar + d2Ubar ) + beta
    dPSI0bottom = - dUbar(1)
    dPSI0top    = - dUbar(NSz)
    
    ! ###################################################################################

    ! matrices for the PV inversion problem
    allocate(A0(NSz,NSz), Ak2(NSz,NSz))
    
    A0(1,:)=T1(1,:)
    do i=2,NSz-1
        A0(i,:)=f0**2/N2bar(i)*( T2(i,:) + (drhobar(i)/rhobar(i) - dN2bar(i)/N2bar(i)) * T1(i,:) )
    end do
    A0(NSz,:)=T1(NSz,:)

    Ak2(1,:)=0
    Ak2(2:NSz-1,:)=-T(2:NSz-1,:)
    Ak2(NSz,:)=0

    ! ###################################################################################

    ! creating the perturbation fields and preparing the simulations to start
    allocate(Q_hat(NSz,NSy,NSx), dPSI_bottom_hat(NSy,NSx), dPSI_top_hat(NSy,NSx), phi(NGx,NGy))
    allocate(Q_hat_old(NSz,NSy,NSx), dPSI_bottom_hat_old(NSy,NSx), dPSI_top_hat_old(NSy,NSx))

    if (kx_initial+ky_initial==0) then ! if the wavenumbers are zero only a random forcing is applied
        call random_number(phi)
        phi=(phi-0.5)*sqrt(12.0)
        call grid_2_spectral_XYplane(phi*Initial_perturbation, dPSI_bottom_hat_old)

    else
        phi(NGx,:)= cos(ky_initial*y/Ly*2*pi) * sqrt(2.0)
        do i=1,NGx
            phi(i,:) = cos(kx_initial*x(i)/Lx*2*pi) * phi(NGx,:)          
        end do
        call grid_2_spectral_XYplane(phi*Initial_perturbation, dPSI_bottom_hat_old)

    end if

    dPSI_top_hat_old = 0

    deallocate(phi)

    dPSI_bottom_hat=dPSI_bottom_hat_old
    dPSI_top_hat=dPSI_top_hat_old

    Q_hat_old = 0
    Q_hat = 0

    tempo=0
    iter=0

    ! ###################################################################################

    tempo1=MPI_wtime()

    call Euler(tempo, Q_hat, dPSI_bottom_hat, dPSI_top_hat, maxU, maxV)
    iter=iter+1
    if (rank==0) print*, iter, tempo/3600, maxU, maxV
    
    do while (tempo<Tfinal)
        call Leapfrog(tempo, Q_hat, dPSI_bottom_hat, dPSI_top_hat, Q_hat_old, dPSI_bottom_hat_old, dPSI_top_hat_old, maxU, maxV)
        iter=iter+1
        if (rank==0) print*, iter, tempo/3600, maxU, maxV
        if (mod(iter,Nsave)==0) then
            call savefields(tempo, Q_hat, dPSI_bottom_hat, dPSI_top_hat, iter)
        end if
    end do
    tempo2=MPI_wtime()
    ! if (rank==0) print*, 'Time for the simulation: ', tempo2-tempo1

    call savefields(tempo, Q_hat, dPSI_bottom_hat, dPSI_top_hat, iter)
    
    ! ###################################################################

    call terminate_FFTW    
    call MPI_Finalize(ierr) 

    contains

    function interp1d(x, y, x0) result(y0)
        integer :: i
        real(8), dimension(:), intent(in) :: x, y
        real(8), intent(in) :: x0
        real(8) :: y0

        ! Find the interval in which x0 lies
        if (x0 < x(1)) then
            y0 = y(1)
        else if (x0 > x(size(x))) then
            y0 = y(size(y))
        else
            do i = 1, size(x)-1
                if (x(i) <= x0 .and. x(i+1) >= x0) then
                    y0 = y(i) + (y(i+1) - y(i)) * (x0 - x(i)) / (x(i+1) - x(i))
                    return
                end if
            end do
        end if
    end function interp1d

    subroutine linear(Q_hat, dPSI_bottom_hat, dPSI_top_hat, NL_Q_hat, NL_dPSI_bottom_hat, NL_dPSI_top_hat, maxUl, maxVl)
        ! linear computations
        integer :: ii, jj, info
        integer, dimension(NSz) :: ipiv

        real(8), dimension(NGz,NGx,NGy) :: U, V
        real(8), intent(out) :: maxUl, maxVl
        real(8), dimension(NGz) :: maxVll

        complex(8), dimension(NSz,NSy,NSx), intent(in)  :: Q_hat
        complex(8), dimension(NSz,NSy,NSx), intent(out) :: NL_Q_hat
        complex(8), dimension(NSz,NSy,NSx)              :: PSI_hat, AUX

        complex(8), dimension(NSy,NSx), intent(in)  :: dPSI_bottom_hat, dPSI_top_hat
        complex(8), dimension(NSy,NSx), intent(out) :: NL_dPSI_bottom_hat, NL_dPSI_top_hat
        complex(8), dimension(NSy,NSx)              :: AUXbottom, AUXtop        
        
        maxUl=0
        ! ###################################################################

        ! inversion of the PV
        PSI_hat(1,:,:)=dPSI_bottom_hat
        PSI_hat(NSz,:,:)=dPSI_top_hat
        do ii=2,NSz-1
            PSI_hat(ii,:,:) = 0
            do jj=1,NSz
                PSI_hat(ii,:,:) = PSI_hat(ii,:,:) + T(ii,jj)*Q_hat(jj,:,:)  
            end do
        end do 

        do jj=1,NSx
            do ii=1,NSy
                call zgesv(NSz, 1, A0+k2(ii,jj)*Ak2, NSz, ipiv, PSI_hat(:,ii,jj), NSz, info)
            end do  
            AUX(:,:,jj) = Iunit * kx(jj) * PSI_hat(:,:,jj) ! this is V_hat          
        end do


        call spectral_2_grid_FFT(AUX,V)
        do ii=1,NGz
            maxVll(ii) = sum(V(ii,:,:)**2)
        end do

        do jj=1,NSx
            AUX(:,:,jj)     = Iunit * kx(jj) * Q_hat(:,:,jj)
            AUXbottom(:,jj) = Iunit * kx(jj) * dPSI_bottom_hat(:,jj) * Ubar(1)
            AUXtop(:,jj)    = Iunit * kx(jj) * dPSI_top_hat(:,jj) * Ubar(NSz)
        end do

        call spectral_2_grid_FFT(AUX, U)
                
        call grid_2_spectral_XYplane( dPSI0bottom * V(1,:,:), NL_dPSI_bottom_hat)
        call grid_2_spectral_XYplane( dPSI0top * V(NGz,:,:),  NL_dPSI_top_hat)
        
        NL_dPSI_bottom_hat = NL_dPSI_bottom_hat + AUXbottom - (N2bar(1)*friction_bottom/f0) * k2 * sum(PSI_hat,1)
        NL_dPSI_top_hat    = NL_dPSI_top_hat    + AUXtop
        
        do ii=1,NGz
            V(ii,:,:) = U(ii,:,:)*Ubar_extended(ii) + Q0_extended(ii)*V(ii,:,:)
        end do
        call grid_2_spectral_FFT( V, NL_Q_hat)  
        
        if (rank==0) then
            call MPI_Reduce(MPI_IN_PLACE, maxVll, NGz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
        else
            call MPI_Reduce(maxVll, maxVll, NGz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
        end if

        maxVl=sqrt(maxval(maxVll)/(NGx*NGy))

    end subroutine linear

    subroutine nonlinear(Q_hat, dPSI_bottom_hat, dPSI_top_hat, NL_Q_hat, NL_dPSI_bottom_hat, NL_dPSI_top_hat, maxUl, maxVl)
        ! nonlinear computations
        integer :: ii, jj, info
        integer, dimension(NSz) :: ipiv

        real(8), dimension(NGz,NGx,NGy) :: U, V, DER
        real(8), dimension(NGx,NGy)     :: DERbottom, DERtop, DERb2, DERt2
        real(8), intent(out) :: maxUl, maxVl
        real(8), dimension(NGz) :: maxUll, maxVll

        complex(8), dimension(NSz,NSy,NSx), intent(in)  :: Q_hat
        complex(8), dimension(NSz,NSy,NSx), intent(out) :: NL_Q_hat
        complex(8), dimension(NSz,NSy,NSx)              :: PSI_hat, AUX

        complex(8), dimension(NSy,NSx), intent(in)  :: dPSI_bottom_hat, dPSI_top_hat
        complex(8), dimension(NSy,NSx), intent(out) :: NL_dPSI_bottom_hat, NL_dPSI_top_hat
        complex(8), dimension(NSy,NSx)              :: AUXbottom, AUXtop        
        
        ! ###################################################################

        ! inversion of the PV
        PSI_hat(1,:,:)=dPSI_bottom_hat
        PSI_hat(NSz,:,:)=dPSI_top_hat
        do ii=2,NSz-1
            PSI_hat(ii,:,:) = 0
            do jj=1,NSz
                PSI_hat(ii,:,:) = PSI_hat(ii,:,:) + T(ii,jj)*Q_hat(jj,:,:)  
            end do
        end do 

        do ii=1,NSy
            do jj=1,NSx
                call zgesv(NSz, 1, A0+k2(ii,jj)*Ak2, NSz, ipiv, PSI_hat(:,ii,jj), NSz, info)
            end do
            AUX(:,ii,:) = - Iunit * ky(ii) * PSI_hat(:,ii,:) ! this is U_hat
        end do

        if (rank==0) then
            AUX(:,1,1)=matmul(Tinv,Ubar)*Nx*Ny ! adding the background flow at the zeroth wavenumber
        end if

        call spectral_2_grid_FFT(AUX,U) ! this is the zonal velocity
        do ii=1,NGz
            maxUll(ii) = sum(U(ii,:,:)**2)
        end do

        do jj=1,NSx
            AUX(:,:,jj)     = Iunit * kx(jj) * Q_hat(:,:,jj)
            AUXbottom(:,jj) = Iunit * kx(jj) * dPSI_bottom_hat(:,jj)
            AUXtop(:,jj)    = Iunit * kx(jj) * dPSI_top_hat(:,jj)
        end do

        call spectral_2_grid_FFT(AUX, DER)
        call spectral_2_grid_XYplane(AUXbottom, DERbottom)
        call spectral_2_grid_XYplane(AUXtop, DERtop)
        
        DERb2 = U(1,:,:) * DERbottom ! first term of the NL term in the bottom boundary
        DERt2 = U(NGz,:,:) * DERtop  ! first term of the NL term in the top boundary
        U     = U * DER              ! second term of the NL term for the PV

        do jj=1,NSx
            AUX(:,:,jj) = Iunit * kx(jj) * PSI_hat(:,:,jj) ! this is V_hat
        end do

        call spectral_2_grid_FFT(AUX,V)
        do ii=1,NGz
            maxVll(ii) = sum(V(ii,:,:)**2)
        end do

        do ii=1,NSy
            AUX(:,ii,:)     = Iunit * ky(ii) * Q_hat(:,ii,:)
            AUXbottom(ii,:) = Iunit * ky(ii) * dPSI_bottom_hat(ii,:)
            AUXtop(ii,:)    = Iunit * ky(ii) * dPSI_top_hat(ii,:)
        end do

        if (rank==0) then
            AUX(:,1,1)      =  matmul(Tinv,Q0)*Nx*Ny ! adding the background flow evolution at the zeroth wavenumber
            AUXbottom(1,1)  =  dPSI0bottom*Nx*Ny
            AUXtop(1,1)     =  dPSI0top*Nx*Ny
        end if
        
        call spectral_2_grid_FFT(AUX,DER)
        call spectral_2_grid_XYplane(AUXbottom,DERbottom)
        call spectral_2_grid_XYplane(AUXtop,DERtop)
        
        call grid_2_spectral_FFT(U + V*DER, NL_Q_hat)
        call grid_2_spectral_XYplane(DERb2 + V(1,:,:)*DERbottom, NL_dPSI_bottom_hat)
        call grid_2_spectral_XYplane(DERt2 + V(NGz,:,:)*DERtop,  NL_dPSI_top_hat)
        
        NL_dPSI_bottom_hat = NL_dPSI_bottom_hat - (N2bar(1)*friction_bottom/f0) * k2 * sum(PSI_hat,1)

        if (rank==0) then
            call MPI_Reduce(MPI_IN_PLACE, maxUll, NGz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(MPI_IN_PLACE, maxVll, NGz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
        else
            call MPI_Reduce(maxUll, maxUll, NGz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            call MPI_Reduce(maxVll, maxVll, NGz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
        end if

        maxUl=sqrt(maxval(maxUll/(NGx*NGy)-Ubar**2))
        maxVl=sqrt(maxval(maxVll)/(NGx*NGy))

    end subroutine nonlinear

    subroutine Euler(tempo, Q_hat, dPSI_bottom_hat, dPSI_top_hat, maxUl, maxVl)
        ! Euler step / implicit in the linear terms
        integer :: ii

        real(8), intent(inout) :: tempo
        real(8), intent(out) :: maxUl, maxVl
        
        complex(8), dimension(NSz,NSy,NSx) :: Q_hat, NL_Q_hat
        complex(8), dimension(NSy,NSx) :: dPSI_bottom_hat, dPSI_top_hat, NL_dPSI_bottom_hat, NL_dPSI_top_hat
        
        tempo = tempo + DT
        if (nonlinear_condition) then
            call nonlinear(Q_hat, dPSI_bottom_hat, dPSI_top_hat, NL_Q_hat, NL_dPSI_bottom_hat, NL_dPSI_top_hat, maxUl, maxVl)
        else
            call linear(Q_hat, dPSI_bottom_hat, dPSI_top_hat, NL_Q_hat, NL_dPSI_bottom_hat, NL_dPSI_top_hat, maxUl, maxVl)
        end if

        dPSI_bottom_hat = (dPSI_bottom_hat - DT * NL_dPSI_bottom_hat)/OPE2
        do ii=1,NSz
            Q_hat(ii,:,:) = (Q_hat(ii,:,:) - DT * NL_Q_hat(ii,:,:))/OPE2
        end do
        dPSI_top_hat = (dPSI_top_hat - DT * NL_dPSI_top_hat)/OPE2

    end subroutine Euler

    subroutine Leapfrog(tempo, Q_hat, dPSI_bottom_hat, dPSI_top_hat, Q_hat_old, dPSI_bottom_hat_old, dPSI_top_hat_old,maxUl,maxVl)
        ! semi-implicit Leapfrog/Crank_Nicolson step
        integer :: ii

        real(8), intent(inout) :: tempo
        real(8), intent(out) :: maxUl, maxVl
        
        complex(8), dimension(NSz,NSy,NSx) :: Q_hat, Q_hat_old, NL_Q_hat
        complex(8), dimension(NSy,NSx) :: dPSI_bottom_hat, dPSI_top_hat, dPSI_bottom_hat_old, dPSI_top_hat_old, AUX, &
                                          NL_dPSI_bottom_hat, NL_dPSI_top_hat
        
        tempo = tempo + DT
        if (nonlinear_condition) then
            call nonlinear(Q_hat, dPSI_bottom_hat, dPSI_top_hat, NL_Q_hat, NL_dPSI_bottom_hat, NL_dPSI_top_hat, maxUl, maxVl)
        else
            call linear(Q_hat, dPSI_bottom_hat, dPSI_top_hat, NL_Q_hat, NL_dPSI_bottom_hat, NL_dPSI_top_hat, maxUl, maxVl)
        end if

        AUX=dPSI_bottom_hat
        dPSI_bottom_hat = (dPSI_bottom_hat_old*OPE1 - 2*DT * NL_dPSI_bottom_hat)/OPE2
        dPSI_bottom_hat_old=AUX !+ 0.01 * (dPSI_bottom_hat - 2*AUX + dPSI_bottom_hat_old)
        do ii=1,NSz
            AUX=Q_hat(ii,:,:)
            Q_hat(ii,:,:) = (Q_hat_old(ii,:,:)*OPE1 - 2*DT * NL_Q_hat(ii,:,:))/OPE2
            Q_hat_old(ii,:,:)=AUX !+ 0.01 * (Q_hat(ii,:,:) - 2*AUX + Q_hat_old(ii,:,:))
        end do
        AUX=dPSI_top_hat
        dPSI_top_hat = (dPSI_top_hat_old*OPE1 - 2*DT * NL_dPSI_top_hat)/OPE2
        dPSI_top_hat_old=AUX !+ 0.01 * (dPSI_top_hat - 2*AUX + dPSI_top_hat_old)
    end subroutine Leapfrog
        
    subroutine savefields(tempo, Q_hat, dPSI_bottom_hat, dPSI_top_hat, iteration)
        ! just saving V and ZETA for now
        integer :: ii, jj, info, iteration
        integer, dimension(NSz) :: ipiv

        real(8), intent(in) :: tempo
        real(8), dimension(NGz,NGx,NGy) :: V      

        complex(8), dimension(NSz,NSy,NSx), intent(in)  :: Q_hat
        complex(8), dimension(NSz,NSy,NSx)              :: PSI_hat, AUX
        complex(8), dimension(NSy,NSx), intent(in)  :: dPSI_bottom_hat, dPSI_top_hat
             
        character(len=100) :: str3, str4

        ! ###################################################################

        ! inversion of the PV to get the streamfunction
        PSI_hat(1,:,:)=dPSI_bottom_hat
        PSI_hat(NSz,:,:)=dPSI_top_hat
        do ii=2,NSz-1
            PSI_hat(ii,:,:) = 0
            do jj=1,NSz
                PSI_hat(ii,:,:) = PSI_hat(ii,:,:) + T(ii,jj)*Q_hat(jj,:,:)  
            end do
        end do 
        do ii=1,NSy
            do jj=1,NSx
                call zgesv(NSz, 1, A0+k2(ii,jj)*Ak2, NSz, ipiv, PSI_hat(:,ii,jj), NSz, info)
            end do
            ! AUX(:,ii,:)=-Iunit*ky(ii)*PSI_hat(:,ii,:) ! this is U_hat
        end do
        ! if (rank==0) then
        !     AUX(:,1,1)=matmul(Tinv,Ubar)*Nx*Ny ! adding the background flow at the zeroth wavenumber
        ! end if

        ! call spectral_2_grid_FFT(AUX,U) ! this is the zonal velocity

        ! ###################################################################

        do jj=1,NSx
            AUX(:,:,jj) = Iunit * kx(jj) * PSI_hat(:,:,jj) ! this is V_hat
        end do

        call spectral_2_grid_FFT(AUX,V)

        write(str3,'(I0)') rank
        write(str4,'(I0)') iteration
        open (unit=10, file=trim('./')//adjustl(trim(fileName_root))//trim('/')//adjustl(trim(fileName_root))//&
        trim('_V_')//trim(str4)//trim('_proc_')//trim(str3)//trim('.bin'), form='unformatted', access='stream', action='write')
        write(10) tempo
        write(10) real(NGz,kind=8), real(NGx,kind=8), real(NGy,kind=8)
        write(10) (1-get_z(NGz))/2*Lz
        write(10) x
        write(10) y
        write(10) V
        close(10)

        ! ###################################################################

        do jj=1,NSz
           PSI_hat(jj,:,:) = - PSI_hat(jj,:,:)  * k2
        end do

        call spectral_2_grid_FFT(PSI_hat,V)

        open (unit=10, file=trim('./')//adjustl(trim(fileName_root))//trim('/')//adjustl(trim(fileName_root))//&
        trim('_ZETA_')//trim(str4)//trim('_proc_')//trim(str3)//trim('.bin'), form='unformatted', access='stream', action='write')
        write(10) tempo
        write(10) real(NGz,kind=8), real(NGx,kind=8), real(NGy,kind=8)
        write(10) (1-get_z(NGz))/2*Lz
        write(10) x
        write(10) y
        write(10) V
        close(10)

    end subroutine savefields

end program QG_plane
