import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Utility function for making plots ------------------------------------------------------------
def fig_phase_space(x,y,nbins,xlabel,ylabel):
    # plt.rcdefaults() # return to default style for fonts in figures
    # plt.rcParams['font.family'] = 'serif'
    # plt.rcParams['font.serif'] = ['Times New Roman']
    # plt.rcParams['font.size'] = 12 
    
    fig = plt.figure(figsize=(6,5))
    gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
                           wspace=0.05, hspace=0.05)

    # Main figure
    ax_main = fig.add_subplot(gs[1, 0])
    h = ax_main.hist2d(x, y, bins=nbins, cmin=1)
    ax_main.set_xlabel(xlabel)
    ax_main.set_ylabel(ylabel)
    ax_main.grid(True)
    
    # Horizontal projection
    ax_x = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_x.hist(x, bins=nbins[0], color='gray')
    #ax_x.axis('off') 
    ax_x.tick_params(labelbottom=False)
    ax_x.grid(True)
    
    # Vertical projection
    ax_y = fig.add_subplot(gs[1, 1], sharey=ax_main)
    ax_y.hist(y, bins=nbins[1], orientation='horizontal', color='gray')
    # ax_y.axis('off') 
    ax_y.tick_params(labelleft=False)
    ax_y.grid(True)

    plt.show()
    return fig, ax_main, ax_x, ax_y
# -----------------------------------------------------------------------------------------------

data = np.loadtxt('./outputs/9Li_8Li-GS_H4C2.txt', delimiter='\t')

xf0, pf0, xf1, pf1, t0, pc_f, Ef = data.T

xf = np.column_stack((xf0, xf1))
pf = np.column_stack((pf0, pf1, np.sqrt(1-pf0**2-pf1**2)))

theta_lab = np.arctan(np.sqrt(pf[:,0]**2+pf[:,1]**2)/pf[:,2])*1e3

mask = Ef > 80
# mask = Ef > 0

nbins = [50,50]
_,_,_,_ = fig_phase_space(theta_lab[mask],Ef[mask], nbins, r"$\theta_{lab}$ (mrad)", "Kinetic Energy (MeV)")
