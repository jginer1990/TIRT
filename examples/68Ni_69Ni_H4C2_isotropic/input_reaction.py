reaction_data = {
    
    'target_composition' : { 
        'layer': [ ['H', 4], ['C', 2] ] , # Composition of the target layer: ['Atom symbol', stoichiometry]
        'density': 0.79867,  # (g/cm^3)
        'thickness': 2e-6  # (m)
    },
       
    'mode' : 'binary',
    
    'projectile' : '68Ni',
    'mass_p' : 67.91655, # optional (amu)
    'ejectile' : '69Ni',
    'mass_e' : 68.92029, # optional (amu)
    'target' : '2H', # (either ion or elementary p, n or e) optional if mass_t is defined
    'mass_t' : 2.01355, # optional (amu) if target_atom is defined
    'recoil' : 'p', # (either ion or elementary p, n or e) optional if mass_r is defined
    'mass_r' : 1.0072764665789, # optional (amu) if recoil is defined
    
    # 'data_XS' : '', # optional (file route with angle_CM(deg)-XS(mb/sr))
    
    # 'max_angle' : 2, # optional (deg): maximum accepted angle of transmitted ions

}