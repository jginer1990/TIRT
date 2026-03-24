beam_config = {
    'mode': 'covariance_6d',
    'N_ions': 100000,
    # 'N_transmitted': 1000, # optional: minimum of transmitted ions needed
    'E0/A' : 10, # (MeV per nucleon)
    # 'E0' : , # optional (MeV)
    # 'pc0': , # optional (MeV)
    'matrix': [
        [0,   0.0,   0.0,    0.0,    0.0,    0.0],
        [0.0,    0,  0.0,    0.0,    0.0,    0.0],
        [0.0,    0.0,   0,    0.0,   0.0,    0.0],
        [0.0,    0.0,   0.0,    0,   0.0,    0.0],
        [0.0,    0.0,   0.0,    0.0,    0,   0.0],        
        [0.0,    0.0,   0.0,    0.0,    0.0,    0]
    ]
}
