reaction_data = {
    
    'target_composition' : { 
        'layer': [ ['H', 4], ['C', 2] ] , # Composition of the target layer: ['Atom symbol', stoichiometry]
        'density': 0.79867,  # (g/cm^3)
        'thickness': 1.25e-6  # (m)
    },
       
    'mode' : 'binary',
    
    'projectile' : '19Ne',
    'mass_p' : 19.001880906, # optional (amu)
    'ejectile' : '20Na',
    'mass_e' : 20.007354301, # optional (amu)
    'target' : '2H', # (either ion or elementary p, n or e) optional if mass_t is defined
    'mass_t' : 2.01355, # optional (amu) if target_atom is defined
    'recoil' : 'n', # (either ion or elementary p, n or e) optional if mass_r is defined
    'mass_r' : 1.00866491606, # optional (amu) if recoil is defined
    
    'data_XS' : '19Ne_20Na_CM.txt', # optional (file route with angle_CM(deg)-XS(mb/sr))
    
    'max_angle' : 2, # optional (deg): maximum accepted angle of transmitted ions

}