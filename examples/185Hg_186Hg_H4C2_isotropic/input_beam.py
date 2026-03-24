beam_config = {
    'mode': 'subspaces',
    'N_ions': 10000,
    # 'N_transmitted': 300, # optional: minimum of transmitted ions needed

    'E0/A' : 10, # (MeV per nucleon)
    # 'E0' : , # optional (MeV)
    # 'pc0': , # optional (MeV)
    
    # X(m)-X'(rad) subspace
    'x': {
        'mode': 'rms',
        'rms_pos': 0.43e-3,  # rms in x or y or time
        'rms_div': 4.25e-3  # rms in x' or y' or dp/p
    },
    
    # Y(m)-Y'(rad) subspace
    'y': {
        'mode': 'rms',
        'rms_pos': 0.43e-3,  # rms in x or y or time
        'rms_div': 4.25e-3  # rms in x' or y' or dp/p
    },
    
    # Longitudinal time(s)-dp/p subspace
    'z': {
        'mode': 'rms',
        'rms_pos': 1e-9/2.355,  # rms in x or y or time
        'rms_div': 0.215/100  # rms in x' or y' or dp/p
    }
}