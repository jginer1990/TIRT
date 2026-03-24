reaction_data = {
    
    'target_composition' : { 
        'layer': [ ['Pb', 1] ] , # Composition of the target layer: ['Atom symbol', stoichiometry]
        'density': 11.34,  # (g/cm^3)
        'thickness': 1.25e-6  # (m)
    },
       
    'mode' : 'binary',
    
    'projectile' : '9Li',
    'mass_p' : 9.02514, # optional (amu)
    'ejectile' : '10Be',
    'mass_e' : 10.01353469, # optional (amu)
    'target' : '208Pb', # (either ion or elementary p, n or e) optional if mass_t is defined
    'mass_t' : 207.97665, # optional (amu) if target_atom is defined
    'recoil' : '207Tl', # (either ion or elementary p, n or e) optional if mass_r is defined
    'mass_r' : 206.97742, # optional (amu) if recoil is defined
    
    # 'data_XS' : '', # optional (file route with angle_CM(deg)-XS(mb/sr))
    
    # 'max_angle' : 2, # optional (deg): maximum accepted angle of transmitted ions

}