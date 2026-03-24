import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import constants
import os, re
from mendeleev import isotope, element # pip install mendeleev

amu = constants.physical_constants['atomic mass constant energy equivalent in MeV'][0]


from utils import define_target, initialize_beam, trim_transport_ions


from input_reaction import reaction_data
from input_beam import beam_config

'''
Main code

'''

print('Define target for SRIM-TRIM')
mytarget, mylayer, material = define_target(reaction_data['target_composition'])
print(material)
print(mylayer)


def ion_identifier(ion):
        match = re.match(r"(\d+)([a-zA-Z]+)", ion)
        if match:
            A = int(match.group(1))
            symbol = match.group(2)
            iso = isotope(symbol, A)
            if iso:
                mass = iso.mass # (amu)
            else:
                mass = 0
            return A, symbol, mass
        else:
            raise ValueError('Wrong format in Ion name (A+Symbol, e.g. "12C")')

# Experiment Name
A_p, ion_p, Mp = ion_identifier(reaction_data['projectile'])
if 'mass_p' in reaction_data:
    Mp = reaction_data['mass_p']
Mp = Mp*amu

experiment_name = f'{A_p}{ion_p}_primary_{material}'

# Maximum number of ions to simulate per batch to avoid memory overload
MAX_BATCH_SIZE = 10000

target_transmitted = beam_config.get('N_transmitted')
target_initial = beam_config.get('N_ions')
max_angle = reaction_data.get('max_angle')
print(f'* Simulation launch summary: {target_initial} ions  --  target at {target_transmitted} transmitted ions ')

# Accumulators for the final transmitted ions
Ef_all = []
xf_all = []
pf_all = []
tf_all = []

total_initialized = 0
total_transmitted = 0

# Get thickness
thickness = reaction_data.get('target_composition')['thickness']

while True:
    # --- 1. Determine batch size (N_batch) ---
    if target_transmitted is not None:
        # Stop condition for N_transmitted mode
        if total_transmitted >= target_transmitted:
            break
            
        remaining_transmitted = target_transmitted - total_transmitted
        
        if total_initialized == 0 or total_transmitted == 0:
            # First iteration or zero transmission so far: use target_initial or MAX_BATCH_SIZE
            N_batch = min([MAX_BATCH_SIZE, target_initial])
        else:
            # Estimate needed initial ions based on current transmission rate
            transmission_rate = total_transmitted / total_initialized
            if transmission_rate == 1:
                estimated_needed = int((remaining_transmitted / transmission_rate))
            else:
                # Add a 5% safety margin to the estimation to minimize extra iterations
                estimated_needed = int((remaining_transmitted / transmission_rate) * 1.05)
            # Ensure N_batch is at least 1 and does not exceed MAX_BATCH_SIZE
            N_batch = max(1, min(MAX_BATCH_SIZE, estimated_needed))
    else:
        # Stop condition for N_ions mode
        if total_initialized >= target_initial:
            break
            
        remaining_initial = target_initial - total_initialized
        N_batch = min(MAX_BATCH_SIZE, remaining_initial)

    print(f'\n--- Starting new batch with N = {N_batch} ions ---')

    # --- 2. Initialize and setup target interactions ---
    x0, p0, t0, E0 = initialize_beam(N_batch, beam_config, Mp, A_p)
        
    # --- 3. SRIM-TRIM Transport of projectile ion up to thickness ###
    print('SRIM-TRIM Transport of projectile ion')

    beam_in_proj = np.column_stack([
        E0,                       # Energy [MeV]
        x0[:, 0:2],               # x, y [m]
        np.zeros_like(E0),        # z [m]
        p0                        # px, py, pz
    ])

    Ef, xf, pf, id2 = trim_transport_ions(
        beam_array=beam_in_proj,
        target=mytarget,
        ion_symbol=ion_p,
        ion_mass=Mp/amu,
        experiment_name=experiment_name
    )
    
    # Update the z-coordinate to the reaction depth
    if len(xf) > 0:
        xf[:, -1] = thickness

        # Filter accepted angle if max_angle exists 
        if max_angle is not None:
            mask_angle = np.sqrt(pf[:, 0]**2 + pf[:, 1]**2) / pf[:, 2] < np.tan(np.radians(max_angle))
            Ef = Ef[mask_angle]
            xf = xf[mask_angle]
            pf = pf[mask_angle]
            tf = t0[mask_angle]
            print(f'Filtering ions in maximum accepted angle: {len(Ef)} out of {len(mask_angle)} ({len(Ef)/len(mask_angle)*100:.3g}%)')
            
        # Append results and update counters 
        if len(Ef) > 0:
            Ef_all.append(Ef)
            xf_all.append(xf)
            pf_all.append(pf)
            tf_all.append(t0[:len(Ef)]) # (this code does not calculate changes in the 'time' coordinate)
            
        transmitted_in_batch = len(Ef)
    else:
        transmitted_in_batch = 0
        print('Warning: No projectile ions reached the reaction depth in this batch.')

    total_initialized += N_batch
    total_transmitted += transmitted_in_batch
    print(f'Batch completed: {transmitted_in_batch} ions transmitted. Total transmitted so far: {total_transmitted}')

# --- 6. Final aggregation and formatting ---
if len(Ef_all) > 0:
    Ef_final = np.concatenate(Ef_all)
    xf_final = np.concatenate(xf_all)
    pf_final = np.concatenate(pf_all)
    tf_final = np.concatenate(tf_all)
    
    # # If using N_transmitted, we might have overshot the exact number needed.
    # # Slice the arrays to match exactly N_transmitted.
    # if target_transmitted is not None and len(Ef_final) > target_transmitted:
    #     Ef_final = Ef_final[:target_transmitted]
    #     xf_final = xf_final[:target_transmitted]
    #     pf_final = pf_final[:target_transmitted]
    #     tf_final = tf_final[:target_transmitted]
else:
    Ef_final, xf_final, pf_final = np.array([]), np.array([]), np.array([])

print(f'\nSimulation fully completed: {len(Ef_final)} ions transmitted out of {total_initialized} initialized.')

Ef, xf, pf, tf = Ef_final, xf_final, pf_final, tf_final


# ------------------------------- Plot results -------------------------------------

# Utility function for making plots
def fig_phase_space(x,y,nbins,xlabel,ylabel,fig_file):
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

    plt.savefig(fig_file)
    return fig, ax_main, ax_x, ax_y

nbins = [50,50]

_,_,_,_ = fig_phase_space( np.arctan(np.sqrt(pf[:,0]**2+pf[:,1]**2)/pf[:,2])*1e3, 
                          Ef, 
                          nbins, 
                          r"$\theta_{lab}$ (mrad)", 
                          "Kinetic Energy (MeV)", 
                          f'./outputs/{experiment_name}_angle_Ekin.png')

# ---------------------------- Export data to file ------------------------------------

pc_f = np.sqrt(Ef*(Ef+2*Mp))
data = np.column_stack((xf[:,0], pf[:,0], xf[:,1], pf[:,1], tf, pc_f, Ef))

header = (
    f"# {A_p}{ion_p} {Mp/amu} u\n"
    f"# {'x(m)':<15} {'px(rad)':<15} {'y(m)':<15} {'py(rad)':<15} {'t(s)':<15} {'p(MeV/c)':<15} {'E(MeV)':<15}"
)
# Save data in txt file
np.savetxt(f'./outputs/{experiment_name}.txt', data, header=header, fmt='%-15.8e', comments='', delimiter='\t')