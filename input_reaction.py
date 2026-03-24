# reaction_data = {
    
#     'target_composition' : { 
#         'layer': [ ['H', 4], ['C', 2] ] , # Composition of the target layer: ['Atom symbol', stoichiometry]
#         'density': 0.79867,  # (g/cm^3)
#         'thickness': 1.25e-6  # (m)
#     },
       
#     'mode' : 'binary',
    
#     'projectile' : '9Li',
#     'mass_p' : 9.02514, # optional (amu)
#     'ejectile' : '10Be',
#     'mass_e' : 10.01134, # optional (amu)
#     'target' : '2H', # (either ion or elementary p, n or e) optional if mass_t is defined
#     'mass_t' : 2.0141, # optional (amu) if target_atom is defined
#     'recoil' : 'n', # (either ion or elementary p, n or e) optional if mass_r is defined
#     'mass_r' : 1.00866491606, # optional (amu) if recoil is defined
    
#     'data_XS' : './Data_Li9/9Li_10Be_CM.txt', # optional (file route with angle_CM(deg)-XS(mb/sr))
    
#     'max_angle' : 5, # optional (deg): maximum accepted angle of transmitted ions

# }

# -------------------------------------------------------------

# reaction_data = {
    
#     'target_composition' : { 
#         'layer': [ ['Pb', 1] ] , # Composition of the target layer: ['Atom symbol', stoichiometry]
#         'density': 11.34,  # (g/cm^3)
#         'thickness': 1.25e-6  # (m)
#     },
       
#     'mode' : 'binary',
    
#     'projectile' : '9Li',
#     'mass_p' : 9.02514, # optional (amu)
#     'ejectile' : '10Be',
#     'mass_e' : 10.01353469, # optional (amu)
#     'target' : '208Pb', # (either ion or elementary p, n or e) optional if mass_t is defined
#     'mass_t' : 207.97665, # optional (amu) if target_atom is defined
#     'recoil' : '207Tl', # (either ion or elementary p, n or e) optional if mass_r is defined
#     'mass_r' : 206.97742, # optional (amu) if recoil is defined
    
#     # 'data_XS' : '', # optional (file route with angle_CM(deg)-XS(mb/sr))
    
#     # 'max_angle' : 2, # optional (deg): maximum accepted angle of transmitted ions

# }

# -------------------------------------------------------------

# reaction_data = {
    
#     'target_composition' : { 
#         'layer': [ ['H', 4], ['C', 2] ] , # Composition of the target layer: ['Atom symbol', stoichiometry]
#         'density': 0.79867,  # (g/cm^3)
#         'thickness': 1.25e-6  # (m)
#     },
       
#     'mode' : 'lab_data',
    
#     'projectile' : '185Hg',
#     # 'mass_p' : , # optional (amu)
#     'ejectile' : '186Hg',
#     # 'mass_e' : , # optional (amu)
    
#     # 'data_XS' : '.txt', # optional (file route with angle_lab(deg)-XS(mb/sr))
#     'gaussfit_XS' : {
#         'A': 1511198.285 , # [mb/sr]
#         'x0': 0 , # centroid [deg]
#         'a0': 0.255 , # sigma width [deg]
#         'range_theta' : [0, 1] # [deg]
#     },
    
#     # 'data_E' : '.txt', # optional (file route with E_lab(MeV)-pdf(a.u.))
#     'gaussfit_E' : {
#         'A': [0.023508,0.022125], # A [a.u.] or normalized
#         'x0': [1821.3473,1843.6866] , # centroid [MeV]
#         'a0': [8.5000,9.0000] , # sigma width [MeV]
#     },
    
#     'Ep_ref': 1850, #  [MeV] Mandatory: reference energy of the projectile ion that corresponds to the provided energy spectrum data
    
#     # 'max_angle' : 2, # optional (deg): maximum accepted angle of transmitted ions

# }

# -------------------------------------------------------------

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