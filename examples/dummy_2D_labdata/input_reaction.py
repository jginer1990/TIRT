reaction_data = {
    
    'target_composition' : { 
        'layer': [ ['H', 4], ['C', 2] ] , # Composition of the target layer: ['Atom symbol', stoichiometry]
        'density': 0.79867,  # (g/cm^3)
        'thickness': 0  # (m)
    },
       
    'mode' : 'lab_data',
    
    'projectile' : '68Ni',
    # 'mass_p' : 9.02514, # optional (amu)
    'ejectile' : '69Ni',
    # 'mass_e' : 10.01134, # optional (amu)
    
    'data_XS_E' : 'dummy_2D_data.txt', # file route with angle_CM(deg)-XS(mb/sr)
    'Ep_ref': 680, #  [MeV] Mandatory: reference energy of the projectile ion that corresponds to the provided energy spectrum data
    
    # 'max_angle' : 5, # optional (deg): maximum accepted angle of transmitted ions

}