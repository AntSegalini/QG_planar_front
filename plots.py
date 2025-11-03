from load_case import read_grid_BIN_file
import numpy as np
from matplotlib.animation import FFMpegWriter
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid



plt.rcParams['pgf.preamble']  =  r'\usepackage{mathptmx}'  # load times roman font
plt.rcParams['font.family']  =  'serif'  # use serif font as default
plt.rcParams['font.size']  =  20  # use serif font as default
plt.rcParams['text.usetex']  =  True  # enable LaTeX rendering globally

# Adjust if needed
rho0 = 1
f0 = 1e-4

CASE_NAME = './CASE_M6_nu30'
# In case of plot failing check the CASE folder to see if the actual time exist

def compute_Q(U, V, B, dx, dy):
    # Compute gradients
    dUdx, dUdy = np.gradient(U, dx, dy, edge_order=2)
    dVdx, dVdy = np.gradient(V, dx, dy, edge_order=2)
    dBdx, dBdy = np.gradient(B, dx, dy, edge_order=2)

    # Q-vector components (standard QG definition)
    Qx = - (dUdx * dBdx + dVdx * dBdy)
    Qy = - (dUdy * dBdx + dVdy * dBdy)
    return Qx, Qy


def plot_B_ZETA_sequence(foldername, colormap = 'turbo', level=0):
    """Animate B and ZETA slices from read_grid_BIN_file.
    Args:
        foldername (str): data path.
        colormap (str): matplotlib colormap.
        level (int): vertical index.
    """
    
    plt.figure(figsize=(16, 12))  # Increase the figure size for a bigger plot window
    for jj in range(0, 3840, 60):
        time, x, y, z, ZETA, U, V, B = read_grid_BIN_file(str(foldername), jj)
        dy = y[1] - y[0]  # assuming uniform spacing
        dx = x[1] - x[0]
        ax = plt.subplot(211)
        B_slice = B[level, :, :]
        B_gradient_y, B_gradient_x = np.gradient(B_slice, dy, dx, edge_order=2)
        grad_magnitude = np.sqrt(B_gradient_x[:, :]**2 + B_gradient_y[:, :]**2)
        plt.contourf(x/1e3, y/1e3, B_slice.T * 288.15 / 9.81, levels=500, cmap=colormap)  # Transpose for z on y-axis
        plt.colorbar()
        
        plt.subplot(212)
        plt.contourf(x, y, ZETA[level, :, :].T, levels=100, cmap='twilight')  # Increased levels for higher resolution
        plt.colorbar()
        # plt.title('Vorticity (Mid-plane)')
        plt.xlabel('x')
        plt.ylabel('z')

        plt.suptitle(rf'Time: {int(time.item() / 3600)} hours / {time.item()/3600/24:.2f} days')
        plt.pause(0.01)
        plt.clf()
        
# plot_B_ZETA_sequence(CASE_NAME)

def B_video(foldername, colormap='RdBu_r'):
    
    """
    Generate an MP4 video animation of the buoyancy field (potential temperature) at the bottom boundary.
    
    Reads binary simulation output files at regular intervals and creates contour plots of the
    buoyancy field with overlaid velocity vectors. Displays time progression in both hours and days.
    Uses a customizable colormap (default: diverging 'RdBu_r' scheme). Domain coordinates are scaled
    from meters to 1000 km for readability. A persistent colorbar is created on the first frame to
    track buoyancy values throughout the animation.
    
    Parameters
    ----------
    foldername : str
        Directory path containing simulation output binary files
    colormap : str, optional
        Matplotlib colormap name (default: 'RdBu_r' for diverging color scheme)
    
    Returns
    -------
    None
        Outputs: B_bottom.mp4 file showing temporal evolution of bottom buoyancy field
    
    Output Files
    ------------
    B_bottom.mp4
        Video file (300 dpi, 5 fps) showing buoyancy and velocity evolution at bottom boundary
    """
    
    fig, ax = plt.subplots(figsize=(10, 8))
    writer = FFMpegWriter(fps=5, metadata=dict(artist='Me'))
    
    cbar = None  # Colorbar
    contour = None  # For storing contourf object

    with writer.saving(fig, "B_bottom.mp4", dpi=300):
        for jj in range(0, 21600, 60):
            time, x, y, _, U, V, ZETA, B = read_grid_BIN_file(str(foldername), jj)

            ax.clear()  # Clear the axes

            # Plot new contour
            contour = ax.contourf(x / 1e6, y / 1e6, B[0, :, :].T* 288.15 / 9.81, levels=500, cmap=colormap)

            # Plot quiver
            ax.quiver(x[::5] / 1e6, y[::5] / 1e6, U[0, ::5, ::5].T, V[0, ::5, ::5].T, color='black', scale=500, alpha=0.5)

            ax.set_title(f'BOTTOM - Time: {int(time.item() / 3600)} hours - {time.item() / 3600 / 24:.1f} days')
            ax.set_xlabel('x (1000 km)')
            ax.set_ylabel('y (1000 km)')
            ax.set_ylim(-2, 2)
            # Only add colorbar for the first frame, but keep it for the others
            if jj == 0:
                cbar = fig.colorbar(contour, ax=ax, orientation='vertical', shrink=0.8, pad=0.02)
                cbar.set_label(r'$\theta$ ($^o$C)')
                cbar.set_ticks(np.linspace(-7, 7, 8))
            # For subsequent frames, do not update or add a new colorbar; keep the first one

            writer.grab_frame()

# B_video(CASE_NAME)

    
def plot_one_field_ZETA(foldername, timestep, colormap = 'twilight', level = 0, projection = 'XY'):
    """
    Plot a single field of ZETA at a specified timestep and projection with pressure contours on it.
    Args:
        foldername (str): data path.
        timestep (int): timestep to plot.
        colormap (str): matplotlib colormap.
        level (int): vertical index for XZ or YZ projection, ignored for XY.
        projection (str): 'XY', 'XZ', or 'YZ'.
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    time, x, y, z, U, V, ZETA, B = read_grid_BIN_file(str(foldername), timestep)
    
    from scipy.integrate import cumulative_trapezoid
    rho0 = 1
    f0 = 0.000119
    # time, x, y, z, U, V, ZETA, B = read_grid_BIN_file(folder, timestep)
    pressure = -rho0 * f0 * cumulative_trapezoid(U, y, axis=2, initial=0)
    
    ax.clear()
    if projection == 'XY':
        contour = ax.contourf(x/1e6, y/1e6, ZETA[level, :, :].T/1.194e-4, levels=500, cmap=colormap)
        ax.quiver(x[::5] / 1e6, y[::5] / 1e6, U[level, ::5, ::5].T, V[level, ::5, ::5].T, color='black', scale=500, alpha=0.5)
        # Add pressure contour lines every 10 units
        pressure_field = pressure[level, :, :].T
        min_p = np.nanmin(pressure_field)
        max_p = np.nanmax(pressure_field)
        levels_p = np.arange(np.floor(min_p / 200) * 200, np.ceil(max_p / 200) * 200 + 200, 200)
        cs = ax.contour(x/1e6, y/1e6, pressure_field, levels=levels_p, colors='red', linewidths=1)
        ax.clabel(cs, inline=True, fontsize=10, fmt='%1.0f')
        ax.set_xlabel('x (1000 km)')
        ax.set_ylabel('y (1000 km)')
        ax.set_title(f'Timestep: {timestep} (Time: {time.item() / 3600 / 24:.2f} days)')
        
    elif projection == 'XZ':
        contour = ax.contourf(x/1e3, z/1e3, ZETA[:, :, level], levels=500, cmap=colormap)
        ax.set_title(f'Timestep: {timestep} (Time: {time.item() / 3600 / 24:.2f} days)')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('z (km)')
    elif projection == 'YZ':
        contour = ax.contourf(y/1e3, z/1e3, ZETA[:, level, :], levels=500, cmap=colormap)
        ax.set_title(f'Timestep: {timestep} (Time: {time.item() / 3600 / 24:.2f} days)')
        ax.set_xlabel('y (km)')
        ax.set_ylabel('z (km)')
    else:
        raise ValueError("Projection must be 'XY', 'XZ', or 'YZ'")

    cbar = fig.colorbar(contour, ax=ax, orientation='vertical', shrink=0.8, pad=0.02)
    cbar.set_label(r'$\zeta$ ($f_0$ units)')
    plt.savefig(f'vort1.png', dpi=300, bbox_inches='tight')
    plt.show()
    
# plot_one_field_ZETA(CASE_NAME, timestep=1440, colormap = 'twilight', level = 0, projection = 'XY')


def time_series_Q_magnitude(foldername, z_level=0):
    """
    Plot time series of the integrated Q magnitude at a specified vertical level.
    """
    times = []
    Q_integrated = []
    max_gradB = []
    for jj in range(0, 15*1440, 60):  # or your timestep range
        time, x, y, z, U, V, ZETA, B = read_grid_BIN_file(foldername, jj)

        dx = x[1] - x[0]
        dy = y[1] - y[0]

        Qx, Qy = compute_Q(U[z_level], V[z_level], B[z_level], dx, dy)
        Qmag = np.sqrt(Qx**2 + Qy**2)

        # Integrate over domain
        Q_total = np.sum(Qmag) * dx * dy

        times.append(time / 86400)  # Convert seconds to days if needed
        Q_integrated.append(np.max(Q_total))  # or np.mean(Q_total) if you prefer
        
        # Compute gradients
        B_gradient_z, B_gradient_x, B_gradient_y = np.gradient(B, z, x, y, edge_order=2)
        
        # Compute gradient magnitude
        gradB_mag = np.sqrt(B_gradient_x**2 + B_gradient_y**2 + B_gradient_z**2)
        
        # Normalize by max value
        gradB_mag= np.max(gradB_mag)  # Or use np.mean(gradB_mag)
        
        # Store
        max_gradB.append(gradB_mag)



    # Normalize Q values by the maximum
    Q_integrated = np.array(Q_integrated)
    Q_integrated_normalized = Q_integrated / np.max(Q_integrated)
    max_gradB = np.array(max_gradB)
    max_grabB_normalized = max_gradB / np.max(max_gradB)
    plt.figure(figsize=(10,6))
    plt.plot(times[::5], Q_integrated_normalized[::5], marker='s', color='black', markersize=4)
    # plt.plot(times, max_grabB_normalized, marker='o', color = 'blue')
    plt.xlabel('Time (days)')
    plt.ylabel(r'$|\mathbf{Q}|$')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('Q_magnitude_time_series.pdf', dpi=500, bbox_inches='tight')
    plt.show()

# time_series_Q_magnitude(CASE_NAME, z_level=0)


def plot_two_field_ZETA(foldername, timestep, colormap='twilight', levels=(0, -1), projection='XY'):
    """ Plot two fields of ZETA at specified levels side by side.
    Args:
        foldername (str): data path.
        timestep (int): timestep to plot.
        colormap (str): matplotlib colormap.
        levels (tuple): vertical indices for the two plots.
        projection (str): 'XY', 'XZ', or 'YZ'.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    time, x, y, z, U, V, ZETA, B = read_grid_BIN_file(str(foldername), timestep)
    
    for i, (ax, level) in enumerate(zip(axes, levels)):
        zeta_normalized = ZETA[level, :, :] / np.max(np.abs(ZETA[level, :, :]))
        contour = ax.contourf(x / 1e6, y / 1e6, zeta_normalized.T, levels=500, cmap=colormap)
        print(f"Level {level}: Max value = {np.max(ZETA[level, :, :]/1.26e-4)}, Min value = {np.min(ZETA[level, :, :]/1.26e-4)}")
        ax.set_xlabel('x (1000 km)')
        # Add text at the left top corner
        # ax.set_ylim(-1.5, 1.5)
        # ax.set_xlim(1.5, 6.5)
        
        if i == 1:
            ax.set_yticklabels([])  # Remove y axis numbers but keep the ticks
            ax.text(0.01, 0.98, r'\textbf{TOP}', transform=ax.transAxes, fontsize=24, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=0, edgecolor='none'))
        else:
            ax.set_ylabel('y (1000 km)')
            ax.text(0.01, 0.98, r'\textbf{BOTTOM}', transform=ax.transAxes, fontsize=24, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=0, edgecolor='none'))
    plt.show()
    
# plot_two_field_ZETA(CASE_NAME, timestep = 0, colormap='twilight_shifted', levels=(0, -1), projection='XY')


def plot_zonalwind(foldername, timestep, colormap1='bwr'):
    """
    Plot zonal wind at a specific timestep.
    """
    time, x, y, z, U, V, ZETA, B = read_grid_BIN_file(str(foldername), timestep)
    
    plt.figure(figsize=(22, 6))
    
    # Plot 1: Zonal wind in XZ plane
    absmax = np.max(np.abs(U))
    vmin, vmax = -absmax, absmax

    plt.subplot(131)
    plt.contourf(y / 1e6, z / 1e3, U[:, U.shape[2] // 2, :], levels=100, cmap=colormap1, vmin=vmin, vmax=vmax)
    plt.ylabel('z (km)')
    plt.xlabel('y (1000 km)')
    plt.colorbar(label=r'$U$ (m/s)', ticks=np.arange(np.ceil(vmin / 10) * 10, vmax + 10, 10))

    # Plot 2: Zonal wind in XY plane (top level)
    plt.subplot(132)
    plt.contourf(x / 1e6, y / 1e6, U[-1, :, :].T, levels=100, cmap=colormap1, vmin=vmin, vmax=vmax)
    plt.ylabel('y (1000 km)')
    plt.xlabel('x (1000 km)')
    plt.colorbar(label=r'$U$ (m/s)', ticks=np.arange(np.ceil(vmin / 10) * 10, vmax + 10, 10))

    # Plot 3: Zonal wind in XY plane (bottom level)
    plt.subplot(133)
    plt.contourf(x / 1e6, y / 1e6, U[0, :, :].T, levels=100, cmap=colormap1, vmin=vmin, vmax=vmax)
    plt.ylabel('y (1000 km)')
    plt.xlabel('x (1000 km)')
    plt.colorbar(label=r'$U$ (m/s)', ticks=np.linspace(vmin, vmax, 43))

    plt.tight_layout()
    plt.show()

# plot_zonalwind(CASE_NAME, 0, colormap1='bwr')


def plot_one_field_ZETA_mult(foldername, timestep=[0, 3, 6, 12,12,12], colormap='twilight', level=0):
    """ Plot multiple fields of ZETA at specified timesteps in a 3x2 grid."""
    fig, axes = plt.subplots(3, 2, figsize=(12, 14), constrained_layout=True)
    
    # First subplot
    time, x, y, z, U, V, ZETA, B = read_grid_BIN_file(str(foldername), timestep[0])
    ax = axes[0, 0]
    zeta_normalized = ZETA[level, :, :] / np.max(np.abs(ZETA[level, :, :]))
    cmap = plt.cm.get_cmap(colormap)
    norm = plt.Normalize(vmin=-1, vmax=1)
    contour = ax.contourf(x / 1e6, y / 1e6, zeta_normalized.T, levels=500, cmap=cmap, norm=norm)
    ax.quiver(x[::7] / 1e6, y[::7] / 1e6, U[level, ::7, ::7].T, V[level, ::7, ::7].T, color='black', scale=500, alpha=0.6)
    print(f"Max value {int(time.item() / 3600 / 24)} days : {np.max(ZETA[level, :, :]/1.194e-4)}, Min value: {np.min(ZETA[level, :, :]/1.194e-4)}")
    ax.text(0.02, 0.98, f'{int(time.item() / 3600 / 24)} days', transform=ax.transAxes, 
            fontsize=25, verticalalignment='top', horizontalalignment='left',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    # ax.set_xlabel('x (km)')
    ax.set_ylim(-1.500, 1.500)
    ax.set_ylabel('y (1000 km)')
    ax.set_xticklabels([])  # Remove numbers from x-axis but keep the ticks
    
    # Second subplot
    time, x, y, z, U, V, ZETA, B = read_grid_BIN_file(str(foldername), timestep[1])
    ax = axes[0, 1]
    zeta_normalized = ZETA[level, :, :] / np.max(np.abs(ZETA[level, :, :]))
    contour = ax.contourf(x / 1e6, y / 1e6, zeta_normalized.T, levels=500, cmap=cmap, norm=norm)
    ax.quiver(x[::7] / 1e6, y[::7] / 1e6, U[level, ::7, ::7].T, V[level, ::7, ::7].T, color='black', scale=500, alpha=0.6)
    print(f"Max value {int(time.item() / 3600 / 24)} days : {np.max(ZETA[level, :, :]/1.194e-4)}, Min value: {np.min(ZETA[level, :, :]/1.194e-4)}")
    ax.text(0.02, 0.98, f'{int(time.item() / 3600 / 24)} days', transform=ax.transAxes, 
            fontsize=25, verticalalignment='top', horizontalalignment='left',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    ax.set_ylim(-1.500, 1.500)

    ax.set_xticklabels([])  
    ax.set_yticklabels([])  
    
    # Third subplot
    time, x, y, z, U, V, ZETA, B = read_grid_BIN_file(str(foldername), timestep[2])
    ax = axes[1, 0]
    zeta_normalized = ZETA[level, :, :] / np.max(np.abs(ZETA[level, :, :]))
    contour = ax.contourf(x / 1e6, y / 1e6, zeta_normalized.T, levels=500, cmap=cmap, norm=norm)
    ax.quiver(x[::7] / 1e6, y[::7] / 1e6, U[level, ::7, ::7].T, V[level, ::7, ::7].T, color='black', scale=500, alpha=0.6)
    print(f"Max value {int(time.item() / 3600 / 24)} days : {np.max(ZETA[level, :, :]/1.194e-4)}, Min value: {np.min(ZETA[level, :, :]/1.194e-4)}")
    ax.text(0.02, 0.98, f'{int(time.item() / 3600 / 24)} days', transform=ax.transAxes, 
            fontsize=25, verticalalignment='top', horizontalalignment='left',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    ax.set_ylim(-1.500, 1.500)
    # ax.set_xlabel('x (km)')
    ax.set_xticklabels([])  
    ax.set_ylabel('y (1000 km)')
    
    # Fourth subplot
    time, x, y, _, U, V, ZETA, B = read_grid_BIN_file(str(foldername), timestep[3])
    ax = axes[1, 1]
    zeta_normalized = ZETA[level, :, :] / np.max(np.abs(ZETA[level, :, :]))
    contour = ax.contourf(x / 1e6, y / 1e6, zeta_normalized.T, levels=500, cmap=cmap, norm=norm)
    quiver = ax.quiver(x[::7] / 1e6, y[::7] / 1e6, U[level, ::7, ::7].T, V[level, ::7, ::7].T, color='black', scale=500, alpha=0.6)
    print(f"Max value {int(time.item() / 3600 / 24)} days : {np.max(ZETA[level, :, :]/1.194e-4)}, Min value: {np.min(ZETA[level, :, :]/1.194e-4)}")
    ax.text(0.02, 0.98, f'{int(time.item() / 3600 / 24)} days', transform=ax.transAxes, 
            fontsize=25, verticalalignment='top', horizontalalignment='left',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    ax.set_ylim(-1.500, 1.500)
    # ax.set_xlabel('x (km)')
    # ax.set_ylabel('y (km)')
    ax.set_yticklabels([])  
    ax.set_xticklabels([])  
    
    # Fourth subplot
    time, x, y, _, U, V, ZETA, B = read_grid_BIN_file(str(foldername), timestep[4])
    ax = axes[2, 0]
    zeta_normalized = ZETA[level, :, :] / np.max(np.abs(ZETA[level, :, :]))
    contour = ax.contourf(x / 1e6, y / 1e6, zeta_normalized.T, levels=500, cmap=cmap, norm=norm)
    quiver = ax.quiver(x[::7] / 1e6, y[::7] / 1e6, U[level, ::7, ::7].T, V[level, ::7, ::7].T, color='black', scale=500, alpha=0.6)
    print(f"Max value {int(time.item() / 3600 / 24)} days : {np.max(ZETA[level, :, :]/1.194e-4)}, Min value: {np.min(ZETA[level, :, :]/1.194e-4)}")
    ax.text(0.02, 0.98, f'{int(time.item() / 3600 / 24)} days', transform=ax.transAxes, 
            fontsize=25, verticalalignment='top', horizontalalignment='left',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    ax.set_ylim(-1.500, 1.500)
    ax.set_xlabel('x (1000 km)')
    # ax.set_ylabel('y (km)')
    ax.set_ylabel('y (1000 km)')
    
    # Fourth subplot
    time, x, y, _, U, V, ZETA, B = read_grid_BIN_file(str(foldername), timestep[5])
    ax = axes[2, 1]
    zeta_normalized = ZETA[level, :, :] / np.max(np.abs(ZETA[level, :, :]))
    contour = ax.contourf(x / 1e6, y / 1e6, zeta_normalized.T, levels=500, cmap=cmap, norm=norm)
    quiver = ax.quiver(x[::7] / 1e6, y[::7] / 1e6, U[level, ::7, ::7].T, V[level, ::7, ::7].T, color='black', scale=500, alpha=0.6)
    print(f"Max value {int(time.item() / 3600 / 24)} days : {np.max(ZETA[level, :, :]/1.194e-4)}, Min value: {np.min(ZETA[level, :, :]/1.194e-4)}")
    ax.text(0.02, 0.98, f'{int(time.item() / 3600 / 24)} days', transform=ax.transAxes, 
            fontsize=25, verticalalignment='top', horizontalalignment='left',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    ax.set_ylim(-1.500, 1.500)
    ax.set_xlabel('x (1000 km)')
    # ax.set_ylabel('y (km)')
    ax.set_yticklabels([])  # Remove numbers from x-axis but keep the ticks

    # Add a legend for the arrows
    ax.quiverkey(quiver, 0.85, -0.17, 30, r'30 m/s', labelpos='E', coordinates='axes', fontproperties={'size': 15}, color='black')
    # # Add a single colorbar for all subplots
    # cbar = fig.colorbar(contour, ax=axes, orientation='vertical', shrink=0.8, pad=0.02)
    # cbar.set_label(r'$\zeta$ (normalized)')
    plt.savefig('output_plot.png', dpi=300, bbox_inches='tight')
    plt.show()

# plot_one_field_ZETA_mult(CASE_NAME, timestep=[2*24*30,3*24*30,5*24*30, 6*24*30, 7*24*30, 8*24*30], colormap='twilight_shifted',level = 0)


def plot_one_field_B_mult(foldername, timestep=[0, 3, 6, 12], colormap='twilight', level=0):
    """ Plot multiple fields of B at specified timesteps in a 2x2 grid with Q-vectors."""
    fig, axes = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True)
    
    def mask_quiver(Qx, Qy, threshold):
        Qmag = np.sqrt(Qx**2 + Qy**2)
        mask = Qmag >= threshold
        Qx_masked = np.where(mask, Qx, np.nan)
        Qy_masked = np.where(mask, Qy, np.nan)
        return Qx_masked, Qy_masked

    Q_threshold = 1e-11  # Set your threshold here

    for idx, ax in enumerate(axes.flat):
        if idx >= len(timestep):
            break
        time, x, y, z, U, V, ZETA, B = read_grid_BIN_file(str(foldername), timestep[idx])
        B_normalized = B[level, :, :]
        contour = ax.contourf(x/1e6, y/1e6, B_normalized[:, :].T * 288.15 / 9.81, levels=500, cmap=colormap)
        Qx, Qy = compute_Q(U[level], V[level], B[level], x[1]-x[0], y[1]-y[0])
        Qx_masked, Qy_masked = mask_quiver(Qx, Qy, Q_threshold)
        # Uncomment below to plot Q-vectors:
        # ax.quiver(x[::8] / 1e6, y[::8] / 1e6, Qx_masked[::8, ::8].T, Qy_masked[::8, ::8].T, color='white', scale=5e-10, alpha=1, linewidth=3, width=0.005)
        # ax.text(0.02, 0.98, f'{time.item() / 3600 / 24:.2f} days', transform=ax.transAxes, 
        #     fontsize=18, verticalalignment='top', horizontalalignment='left',
        #     bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        # ax.set_ylim(-1.5, 1.5)
        if idx % 2 == 0:
            ax.set_ylabel('y (1000 km)')
        else:
            ax.set_yticklabels([])
        if idx >= 2:
            ax.set_xlabel('x (1000 km)')
        else:
            ax.set_xticklabels([])

    cbar = fig.colorbar(contour, ax=axes, orientation='horizontal', shrink=0.8, pad=0.02, aspect=30, extend='both')
    cbar.set_label(r'$\theta^\prime$ ($^o$C)')
    plt.savefig('output_plot.png', dpi=300, bbox_inches='tight')
    plt.show()

# plot_one_field_B_mult(CASE_NAME, timestep=[2*24*30,3*24*30,5*24*30, 6*24*30], colormap='RdBu_r', level=0)

def plot_sequence_B(foldername, timesteps, colormap='twilight', projection='XY', level=0):
    """
    Plot a 2x5 grid (10 subplots) of B field for 10 different timesteps.
    """
    
    fig, axes = plt.subplots(2, 5, figsize=(14, 7), constrained_layout=True)
    for idx, ax in enumerate(axes.flat):
        if idx >= len(timesteps):
            ax.axis('off')
            continue
        time, x, y, z, U, V, ZETA, B = read_grid_BIN_file(str(foldername), int(timesteps[idx]))

        pressure = -rho0 * f0 * cumulative_trapezoid(U, y, axis=2, initial=0)
        def mask_quiver(Qx, Qy, threshold):
            Qmag = np.sqrt(Qx**2 + Qy**2)
            mask = Qmag >= threshold
            Qx_masked = np.where(mask, Qx, np.nan)
            Qy_masked = np.where(mask, Qy, np.nan)
            return Qx_masked, Qy_masked

        Q_threshold = 1e-30
        U_masked, V_masked = mask_quiver(U[level], V[level], Q_threshold)

        if projection == 'XY':
            contour = ax.contourf(x/1e6, y/1e6, B[level, :, :].T * 288.15 / 9.81, levels=100, cmap=colormap)
            # ax.quiver(x[::8] / 1e6, y[::8] / 1e6, U_masked[::8, ::8].T, V_masked[::8, ::8].T, color='black', scale=5e2, alpha=1, linewidth=1, width=0.004)
            ax.text(0.03, 0.96, fr'{time.item() / 3600 / 24 :.0f} days', transform=ax.transAxes, 
                fontsize=15, verticalalignment='top', horizontalalignment='left',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
            ax.set_xlim(x[0]/1e6, x[-1]/1e6)
            ax.set_ylim(y[0]/1e6, y[-1]/1e6)
            ax.set_aspect('equal')
            ax.set_yticks([-2,-1,0, 1, 2])
            ax.set_xticks([0, 1, 2, 3, 4])
            if idx % 5 == 0:
                ax.set_ylabel('y (1000 km)')
            else:
                ax.set_yticklabels([])
            if idx // 5 == 1:
                ax.set_xlabel('x (1000 km)')
            else:
                ax.set_xticklabels([])
        else:
            ax.axis('off')

    # Add a single colorbar for all subplots
    cbar = fig.colorbar(contour, ax=axes, orientation='vertical', shrink=0.8, pad=0.02)
    cbar.set_label(r'$\theta$ (K)')
    cbar.set_ticks([-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10])
    plt.savefig('sequence_B_10subplots.pdf', dpi=500, bbox_inches='tight')
    plt.show()

# timesteps =  np.array([0 ,1440 *1, 1440 *2, 1440 *4, 1440 *6, 1440 *8, 1440 *10,  1440 *12, 1440 *14, 1440 *9])
# plot_sequence_B(CASE_NAME, timesteps, colormap='RdBu_r',level = 0)



































