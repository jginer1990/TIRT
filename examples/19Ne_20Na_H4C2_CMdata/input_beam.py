beam_config = {
    'mode': 'subspaces',
    'N_ions': 10000,
    'N_transmitted': 10000, # optional: minimum of transmitted ions needed

    'E0/A' : 10, # (MeV per nucleon)
    # 'E0' : , # optional (MeV)
    # 'pc0': , # optional (MeV)
    
    # X(m)-X'(rad) subspace
    'x': {
        'mode': 'twiss',
        'emit': 0.6e-6,  # Emittance (geometric)
        'beta': 1.54,    # Courant-Snyder (Twiss) functions
        'alpha': 0.29
    },
    
    # Y(m)-Y'(rad) subspace
    'y': {
        'mode': 'twiss',
        'emit': 0.6e-6,  # Emittance (geometric)
        'beta': 1.21,    # Courant-Snyder (Twiss) functions
        'alpha': 0.002
    },
    
    # Longitudinal time(s)-dp/p subspace
    'z': {
        'mode': 'rms',
        'rms_pos': 1e-9/2.355,  # rms in x or y or time
        'rms_div': 0.005  # rms in x' or y' or dp/p
    }
}
