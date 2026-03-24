reaction_data = {
    
    'target_composition' : { 
        'layer': [ ['H', 4], ['C', 2] ] , # Composition of the target layer: ['Atom symbol', stoichiometry]
        'density': 0.79867,  # (g/cm^3)
        'thickness': 1.25e-6  # (m)
    },
       
    'mode' : 'lab_data_uncorrelated',
    
    'projectile' : '185Hg',
    # 'mass_p' : , # optional (amu)
    'ejectile' : '186Hg',
    # 'mass_e' : , # optional (amu)
    
    # 'data_XS' : '.txt', # optional (file route with angle_lab(deg)-XS(mb/sr))
    'gaussfit_XS' : {
        'A': 1511198.285 , # [mb/sr]
        'x0': 0 , # centroid [deg]
        'a0': 0.255 , # sigma width [deg]
        'range_theta' : [0, 1] # [deg]
    },
    
    # 'data_E' : '.txt', # optional (file route with E_lab(MeV)-pdf(a.u.))
    'gaussfit_E' : {
        'A': [0.023508,0.022125], # A [a.u.] or normalized
        'x0': [1821.3473,1843.6866] , # centroid [MeV]
        'a0': [8.5000,9.0000] , # sigma width [MeV]
    },
    
    'Ep_ref': 1850, #  [MeV] Mandatory: reference energy of the projectile ion that corresponds to the provided energy spectrum data
    
    # 'max_angle' : 2, # optional (deg): maximum accepted angle of transmitted ions

}