SRIM_DIR = r"C:\Users\jorge\Documents\SRIM-2013"

import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d, RegularGridInterpolator
from scipy import integrate, constants
from scipy.optimize import minimize_scalar, minimize

from mendeleev import isotope, element # pip install mendeleev

import os, re

from srim import TRIM, Ion, Layer, Target

amu = constants.physical_constants['atomic mass constant energy equivalent in MeV'][0]


""" Functions of Utility for TIRT code """

################   REJECTION SAMPLING    #####################################################

def truncated_cauchy_pdf(x, x0, gamma, xmin, xmax):
    """Truncated Cauchy-Lorentz Probability Density Function."""
    # Standard Cauchy CDF: F(x) = (1/pi) * arctan((x-x0)/gamma) + 0.5
    f_cdf = lambda val: np.arctan((val - x0) / gamma) / np.pi + 0.5
    
    # Standard Cauchy PDF
    pdf = 1 / (np.pi * gamma * (1 + ((x - x0) / gamma)**2))
    
    # Normalization factor (Area between xmin and xmax)
    normalization = f_cdf(xmax) - f_cdf(xmin)
    
    # Return PDF if within bounds, else 0
    return np.where((x >= xmin) & (x <= xmax), pdf / normalization, 0)

def truncated_cauchy_rvs(x0, gamma, xmin, xmax, size=1):
    """Random Variates from Truncated Cauchy using Inverse Transform Sampling."""
    f_cdf = lambda val: np.arctan((val - x0) / gamma) / np.pi + 0.5
    
    # Sample uniformly between the CDF values of the boundaries
    low, high = f_cdf(xmin), f_cdf(xmax)
    u = np.random.uniform(low, high, size)
    
    # Inverse CDF: x = x0 + gamma * tan(pi * (u - 0.5))
    return x0 + gamma * np.tan(np.pi * (u - 0.5))


def optimize_cauchy_envelope(interpolator, theta_min, theta_max, total_area_f):
    """
    This function finds the best Cauchy parameters (theta_0, gamma) 
    and the scaling factor M for your specific interpolator and range.
    """
    # 1. Find the peak position (theta_0) within the range to center the Cauchy
    # We sample the interpolator to find the starting maximum
    test_angles = np.linspace(theta_min, theta_max, 1000)
    test_values = interpolator(test_angles)
    theta_0 = test_angles[np.argmax(test_values)]
    
    # # Calculate the total area (integral) of your function for efficiency calculation
    # total_area_f2, _ = quad(interpolator, theta_min, theta_max, limit=100)
    # print(total_area_f, total_area_f2)

    def calculate_M(gamma):
        if gamma <= 0: return np.inf
        
        # Check ratios f/g at multiple points to find the global M
        check_points = np.linspace(theta_min, theta_max, 901)
        g_theta = truncated_cauchy_pdf(check_points, theta_0, gamma, theta_min, theta_max)
        f_theta = interpolator(check_points)
        
        return np.max(f_theta / g_theta)

    # 2. Optimize gamma (width) to minimize the envelope area M
    # Range for gamma search: from 0.001 degrees to twice the total range
    res = minimize_scalar(calculate_M, bounds=(0.001, (theta_max - theta_min) * 2), method='bounded')

    gamma_opt = res.x
    M_opt = res.fun
    efficiency = total_area_f / M_opt # Theoretical acceptance rate
    
    return theta_0, gamma_opt, M_opt, efficiency


def rejection_sampling_cauchy(interpolator, theta_min, theta_max, theta_0, gamma, M):
    """ Rejection sampling algorithm, improved acceptance rate by Cauchy-Lorentz generator """
    while True:
        # 1. Sample theta from the Cauchy distribution using the Inversion Method
        theta_gen = truncated_cauchy_rvs(theta_0, gamma, theta_min, theta_max)
        
        # 2. Evaluate the Cauchy PDF at the generated point
        g_val = truncated_cauchy_pdf(theta_gen, theta_0, gamma, theta_min, theta_max)
        
        # 3. Evaluate your function via the interpolator
        f_val = interpolator(theta_gen)
        
        # 4. Rejection Criterion
        if np.random.random() <= f_val / (M * g_val):
            return theta_gen # Returns the accepted angle in degrees


## ------------- Rejection sampling (optimized) for 2D distribution --------------------- ##

def calculate_2d_volume(theta_grid, energy_grid, prob_matrix):
    """Calculates the total volume (integral) of the 2D probability surface."""
    # Composite trapezoidal rule over energy then theta
    integral_over_E = np.trapz(prob_matrix, x = energy_grid, axis = 1)
    total_volume = np.trapz(integral_over_E, x = theta_grid)
    return total_volume

def create_2d_interpolator_from_file(filename):
    """
    Reads a 3-column file (theta in deg, E in MeV, prob), applies the sin(theta) 
    correction, and creates a RegularGridInterpolator.
    Returns the interpolator, bounds, sampler volume, and total cross-section.
    """
    data = np.loadtxt(filename)
    t_col, e_col, p_col = data[:, 0], data[:, 1], data[:, 2]
    
    # Apply solid angle correction: P' = P * sin(theta)
    # np.sin requires radians, so we convert theta first
    p_col = p_col * np.sin(np.radians(t_col))
    
    # Extract unique, sorted axes to form the regular grid
    u_theta = np.sort(np.unique(t_col))
    u_energy = np.sort(np.unique(e_col))
    
    # Map points to a 2D grid matrix
    prob_matrix = np.zeros((len(u_theta), len(u_energy)))
    for t, e, p in zip(t_col, e_col, p_col):
        t_idx = np.where(u_theta == t)[0][0]
        e_idx = np.where(u_energy == e)[0][0]
        prob_matrix[t_idx, e_idx] = p
        
    interp_2d = RegularGridInterpolator(
        (u_theta, u_energy), 
        prob_matrix, 
        method='linear', 
        bounds_error=False, 
        fill_value=0.0
    )
    
    theta_bounds = (u_theta[0], u_theta[-1])
    energy_bounds = (u_energy[0], u_energy[-1])
    
    # 1. Volume for the Monte Carlo sampler (integrated over degrees)
    sampler_volume = calculate_2d_volume(u_theta, u_energy, prob_matrix)
    
    # 2. Volume for Physics / Total Cross Section (integrated over radians * 2pi)
    physical_volume = calculate_2d_volume(np.radians(u_theta), u_energy, prob_matrix)
    total_cross_section = physical_volume * 2 * np.pi
    
    return interp_2d, theta_bounds, energy_bounds, sampler_volume, total_cross_section

def truncated_cauchy_2d_pdf(x, y, x0, y0, gx, gy, xb, yb):
    """Product of two 1D Truncated Cauchy PDFs for the 2D envelope."""
    return truncated_cauchy_pdf(x, x0, gx, xb[0], xb[1]) * \
           truncated_cauchy_pdf(y, y0, gy, yb[0], yb[1])

def truncated_cauchy_2d_rvs(x0, y0, gx, gy, xb, yb, size=1):
    """Generates independent random variates for x (theta) and y (Energy)."""
    x = truncated_cauchy_rvs(x0, gx, xb[0], xb[1], size)
    y = truncated_cauchy_rvs(y0, gy, yb[0], yb[1], size)
    return x, y

def optimize_cauchy_2d_envelope(interp_2d, x_bounds, y_bounds, total_vol):
    """Finds optimal (theta0, E0, gamma_theta, gamma_E) and scaling factor M."""
    # Find peak on a coarse grid for centering
    tx = np.linspace(x_bounds[0], x_bounds[1], 400)
    ty = np.linspace(y_bounds[0], y_bounds[1], 400)
    X, Y = np.meshgrid(tx, ty)
    Z = interp_2d(np.array([X, Y]).T).T # Ensure RGI evaluation on meshgrid
    
    max_idx = np.unravel_index(np.argmax(Z), Z.shape)
    x0, y0 = tx[max_idx[1]], ty[max_idx[0]]

    def find_M(gammas):
        gx, gy = gammas
        if gx <= 0 or gy <= 0: return np.inf
        g_val = truncated_cauchy_2d_pdf(X, Y, x0, y0, gx, gy, x_bounds, y_bounds)
        return np.max(Z / g_val)

    res = minimize(find_M, [x_bounds[1]/10, y_bounds[1]/10], 
                   bounds=[(0.001, None), (0.001, None)], method='L-BFGS-B')
    
    return (x0, y0), res.x, res.fun, total_vol / res.fun

def rejection_sampling_2d_cauchy(interp_2d, xb, yb, x0, y0, gx, gy, M):
    """Performs Rejection Sampling using the 2D Cauchy envelope."""
    while True:
        xt, yt = truncated_cauchy_2d_rvs(x0, y0, gx, gy, xb, yb)
        # Handle potential array return from rvs
        xt_s, yt_s = (xt[0], yt[0]) if isinstance(xt, np.ndarray) else (xt, yt)
        
        g_val = truncated_cauchy_2d_pdf(xt_s, yt_s, x0, y0, gx, gy, xb, yb)
        f_val = interp_2d([xt_s, yt_s])[0] # RGI evaluation
        
        if np.random.random() <= f_val / (M * g_val):
            return xt_s, yt_s



################   NUCLEAR INTERACTION   #####################################################

def process_reaction(En, pn, reaction_config):
    """
    Calculates ejectile directions (p2) and energy (Ee) from the 
    random angle (theta) generated based on cross-section data
    """
    reaction_type = reaction_config['mode']

    theta, Ee, p2 = [], [], []
    for Ei, p_i in zip(En,pn):
        if reaction_config['mode'] == 'binary':
            xs_interpolator = reaction_config['interpolator']
            range_theta = reaction_config['range_theta']
            xs_max = reaction_config['XS_max']
            cauchy_param = reaction_config['cauchy_param']
            
            theta_CM = rejection_sampling_cauchy(xs_interpolator, 
                                                 range_theta[0], range_theta[1],
                                                 cauchy_param[0], cauchy_param[1],
                                                 xs_max
                                                )
            
            Mp, Me, mt, mr = reaction_config['masses']
            
            theta_, Ee_ = CM_to_LAB(theta_CM, Ei, Mp,mt,Me,mr)  # lab angle [deg] and lab energy [MeV]


        
        elif reaction_type == "lab_data_uncorrelated":
            xs_interpolator = reaction_config['interpolator'][0]
            range_theta = reaction_config['range_theta']
            xs_max = reaction_config['XS_max']
            cauchy_param_xs = reaction_config['cauchy_param_XS']

            theta_ = rejection_sampling_cauchy(xs_interpolator, 
                                               range_theta[0], range_theta[1],
                                               cauchy_param_xs[0], cauchy_param_xs[1],
                                               xs_max
                                               )

            fE_interpolator = reaction_config['interpolator'][1]
            range_E = reaction_config['range_E']
            fE_max = reaction_config['fE_max']
            cauchy_param_E = reaction_config['cauchy_param_E']
            Ep_ref = reaction_config['Ep_ref']
            
            Ee_ = rejection_sampling_cauchy(fE_interpolator, 
                                            range_E[0], range_E[1],
                                            cauchy_param_E[0], cauchy_param_E[1],
                                            fE_max
                                           )
            Ee_ += Ei - Ep_ref # shift with respect to the actual ion energy
            
        elif reaction_type == "lab_data":
            # 1. Extract 2D parameters from config
            interp_2d = reaction_config['interpolator_2d']
            r_theta = reaction_config['range_theta']
            r_E = reaction_config['range_E']
            M_max = reaction_config['M_max']
            centers, gammas = reaction_config['cauchy_param_2d']
            
            # 2. Sample theta and Energy simultaneously
            theta_, Ee_ = rejection_sampling_2d_cauchy(
                interp_2d, r_theta, r_E, 
                centers[0], centers[1], gammas[0], gammas[1], M_max
            )
            
            # 3. Correct energy relative to beam energy
            Ee_ += (Ei - reaction_config['Ep_ref'])
        
        else:
            raise ValueError("Error in reaction type")
            
        p2_ = new_direction(p_i,theta_)
    
        theta.append(theta_)
        Ee.append(Ee_)
        p2.append(p2_)
                
    return np.array(p2), np.array(Ee), np.array(theta)

# --------------------

def new_direction(p0, theta):
    # p0: 3D vector, original direction
    # theta: angle in degrees
        
    rand_values = np.random.rand(2)
    
    # isotropic vector: q(sintheta*sinphi,sintheta*cosphi,costheta)
    qcostheta = 2 * rand_values[0] - 1
    qsintheta = np.sqrt(1 - qcostheta**2)
    qphi      = 2 * np.pi * rand_values[1]

    # new vector u, perpendicular to p0: [u = p0 x q (vector product)]
    u = np.array([
        p0[1] * qcostheta - p0[2] * qsintheta * np.sin(qphi),
        p0[2] * qsintheta * np.cos(qphi) - p0[0] * qcostheta,
        p0[0] * qsintheta * np.sin(qphi) - p0[1] * qsintheta * np.cos(qphi)
    ])
    u = u / np.linalg.norm(u)  # Normalized

    # New direction at an angle 'theta' with respect to p0
    return np.cos(np.radians(theta)) * p0 + np.sin(np.radians(theta)) * u

# ----------------------

def CM_to_LAB(theta_CM, Ep, Mp,mt,Me,mr):
    """
    theta_CM: ejectile angle in CM frame [degrees]
    Ep: projectile kinetic energy [MeV]
    Mp,mt,Me,mr: masses in [MeV/c²] of projectile, target, ejectile and recoil
    """
    
    ECM = mt/(Mp+mt)*Ep
    Q = Mp+mt-Me-mr
    
    gamma = np.sqrt(Mp*Me/mt/mr*(Me+mr)/(Mp+mt)*ECM/(ECM+Q))

    theta = np.degrees( np.arctan2( np.sin(np.radians(theta_CM)), (np.cos(np.radians(theta_CM))+gamma) ) ) 
    # theta: only values between 0 and 180 deg, backscattering if >90 deg (this happens only when "cos(thetaCM)+gamma<0" )
    Ee = mr*(Q+ECM)/(Me+mr)*(1+gamma**2+2*gamma*np.cos(np.radians(theta_CM)))

    return theta, Ee





################   TRIM TRANSPORT   #####################################################



def trim_transport_ions(
    beam_array, # 6 columns
    target, # TRIM target class
    ion_symbol, # string
    ion_mass, # in atomic mass units
    experiment_name="Exp1"
):
    """
    Ion transport with TRIM-pysrim (calculation=4).

    beam_array columns:
    [E(MeV), x(m), y(m), z(m), cosX, cosY, cosZ]

    Returns:
    energies = E(MeV)
    positions = x(m), y(m), z(m)
    directions = cosX, cosY, cosZ
        only for transmitted ions.
    """

    if beam_array.shape[1] != 7:
        raise ValueError("beam_array shape must be (N,7)")

    N = beam_array.shape[0]

    beam = beam_array.copy() # Make a copy to convert units and swap components

    beam[:, 0] *= 1e6        # MeV → eV

    # swap Z axes → X axes from TRIM
    beam[:, 1:4] = beam[:, [3, 2, 1]] * 1e10    # z,y,x [m → Å]
    
    dirs = beam[:, [6, 5, 4]]    # cosZ, cosY, cosX
    norms = np.linalg.norm(dirs, axis=1) # normalize directions
    dirs /= norms[:, None]
    beam[:, 4:7] = dirs

    Z = element(ion_symbol).atomic_number # Get atomic number

    # Write TRIM.DAT file
    print('Write TRIM.DAT file')
    
    trim_dat_path = os.path.join(SRIM_DIR, "TRIM.DAT")

    header_lines = [
        "% TRIM with various Incident Ion Energies/Angles and Depths %",
        "% Top 10 lines are user comments, with line #8 describing experiment. %",
        "% Line #8 will be written into all TRIM output files. %",
        "% Data Table line consist of: EventName(5 char.) + 8 numbers separated by spaces. %",
        "% Event Name: up to 5 characters identifying the particle. %",
        "% Generated automatically by Python interface. %",
        "% cos(X)=1 normal incidence; cos(X)=-1 towards surface %",
        f"% {experiment_name} %",
        "% Event\tAtom\tEnergy\tDepth\t-Lateral-Pos-\t--- Atom Direction ---",
        "% Name\tNumb\t(eV)\tX (A)\tY (A)\tZ (A)\tCos(X)\tCos(Y)\tCos(Z)"
    ] # Header template for TRIM.DAT

    with open(trim_dat_path, "w") as f:

        for line in header_lines:  # write header
            f.write(line + "\n")

        for i in range(N): # write data
            event_name = str(i)[:5]  # ion identifier (max 5 chars) - UNUSED
            row = beam[i]
            f.write(
                f"{event_name:<5}\t"
                f"{Z:>3d}\t"
                f"{row[0]:.6e}\t"
                f"{row[1]:.6e}\t"
                f"{row[2]:.6e}\t"
                f"{row[3]:.6e}\t"
                f"{row[4]:.6e}\t"
                f"{row[5]:.6e}\t"
                f"{row[6]:.6e}\n"
            )

    # Define Ion (max energy needed by TRIM for building tables)
    ion = Ion(ion_symbol, mass=ion_mass, energy=1.05*np.max(beam[:, 0]))

    trim = TRIM(
        target,
        ion,
        number_ions=N,
        calculation=4,
        transmit=1
    )

    trim.run(SRIM_DIR)

    # Read output file TRANSMIT.txt
    print('Read output file TRANSMIT.txt')

    transmit_path = os.path.join(
        SRIM_DIR,
        "SRIM Outputs",
        "TRANSMIT.txt"
    )

    if not os.path.exists(transmit_path):
        raise FileNotFoundError("TRANSMIT.txt NOT FOUND")

    Ek = []
    x, y, z = [], [], []
    px, py, pz = [], [], []
    ids = []
    found_T = False
    with open(transmit_path, 'r') as f:
        for line in f:
            if line.startswith('T'):
                line = line[1:]
                parts = line.split()

                ids.append( int(parts[0]) ) # Ion ID number
                
                Ek.append( float(parts[2])*1e-6 ) # Kinetic energy [MeV]
                x.append( float(parts[5])*1e-10 ) # Lateral X [m] (swap X-Z undone)
                y.append( float(parts[4])*1e-10 ) # Lateral Y [m]
                z.append( float(parts[3])*1e-10 ) # Depth Z [m] (swap X-Z undone)
                px.append( float(parts[8]) ) # Lateral X [m] (swap X-Z undone)
                py.append( float(parts[7]) ) # Lateral Y [m]
                pz.append( float(parts[6]) ) # Depth Z [m] (swap X-Z undone)
                found_T = True

    if not found_T:
        raise ValueError("No ions were transmitted.")

    energies = np.array(Ek)
    positions = np.column_stack([x,y,z])
    directions = np.column_stack([px,py,pz])
    id_ion = np.array(ids)

    return energies, positions, directions, id_ion


# ---------------------------------


def define_target(target_info):
    from atomdata import SRIMdata

    elements_dict = {}

    material = ''
    
    for symbol, stoich in target_info['layer']:
        if symbol in SRIMdata:
            # Combine stoichiometry with SRIMdata
            elements_dict[symbol] = {
                'stoich': stoich,
                'E_d': SRIMdata[symbol]['E_d'],
                'lattice': SRIMdata[symbol]['lattice'],
                'surface': SRIMdata[symbol]['surface']
            }
            material += f"{symbol}{int(stoich)}"
        else:
            print(f"Error: Data for {symbol} not found")
    
    thickness_A = target_info['thickness'] * 1e10 # (Å)
    
    layer = Layer(
        elements_dict, 
        density=target_info['density'], 
        width=thickness_A
    )

    mytarget = Target( [ layer ] ) # pysrim (TRIM) class for target definition
    
    return mytarget, layer, material


################   DATA CONFIGURATION   #####################################################

def gaussfit(x,A,x0,a0):
    # x: angle in [deg] or Energy in [MeV]
    # A,x0,a0: differential cross-section parameters for a "sum-of-gaussians" [a.u.], [deg or MeV], [deg or MeV]

    x = np.atleast_1d(x)  # make sure it is a 1D array
    A = np.atleast_1d(A)
    x0 = np.atleast_1d(x0)
    a0 = np.atleast_1d(a0)
    
    f = A[:, None] * np.exp(-((x[None, :] - x0[:, None])**2) / (2 * a0[:, None]**2))
    
    return np.sum(f, axis=0) # sum over gaussian parameters, return x-like dimension array


# --------------------------

def process_reaction_data(reaction_data):
    
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
    
    
    def element_mass(ele):
        match ele:
            case 'p':
                return constants.physical_constants['proton mass in u'][0]
            case 'n':
                return constants.physical_constants['neutron mass in u'][0]
            case 'e':
                return constants.physical_constants['electron mass in u'][0]
                
    
    
    if reaction_data['mode'] == 'binary':
        
        A_p, ion_p, Mp = ion_identifier(reaction_data['projectile'])
        if 'mass_p' in reaction_data:
            Mp = reaction_data['mass_p']
            
        A_e, ion_e, Me = ion_identifier(reaction_data['ejectile'])
        if 'mass_e' in reaction_data:
            Me = reaction_data['mass_e']
    
        if 'mass_t' in reaction_data:
            mt = reaction_data['mass_t']
        elif 'target' in reaction_data:
            if reaction_data['target'] in ['p', 'n', 'e']:
                mt = element_mass(reaction_data['target'])
            else:
                _, _, mt = ion_identifier(reaction_data['target'])
        else:
            raise ValueError('No target ion or particle is defined')
    
        if 'mass_r' in reaction_data:
            mr = reaction_data['mass_r']
        elif 'target' in reaction_data:
            if reaction_data['recoil'] in ['p', 'n', 'e']:
                mr = element_mass(reaction_data['recoil'])
            else:
                _, _, mr = ion_identifier(reaction_data['recoil'])
        else:
            raise ValueError('No recoil ion or particle is defined')
    
        Mp = Mp*amu 
        mt = mt*amu
        Me = Me*amu
        mr = mr*amu
        
        if 'data_XS' in reaction_data:
            tCM, xsCM = np.loadtxt(reaction_data['data_XS'], unpack=True)
            i_sort = np.argsort(tCM)
            tCM = tCM[i_sort]
            xsCM = xsCM[i_sort]
        else:
            tCM = np.linspace(0, 180, 1801)
            xsCM = np.ones_like(tCM)
        fxsCM = interp1d(tCM, xsCM*np.sin(np.radians(tCM)), kind='linear')
        
        # Check total cross-section
        xx = np.linspace(np.min(tCM),np.max(tCM),2001) # sampling for integral
        xs_tot = integrate.trapezoid(fxsCM(xx), x=np.radians(xx) )*2*np.pi
        print(f"Calculate Total cross section: {xs_tot:.4g} mb")
    
        area_xs = integrate.trapezoid(fxsCM(xx), x=xx )
        theta0_cauchy, gamma_cauchy, MxsCM, eta_rejection = optimize_cauchy_envelope(fxsCM, np.min(tCM), np.max(tCM), area_xs) 
        print('Optimize rejection sampling efficiency')
        print(f'Cauchy parameters: theta_0 = {theta0_cauchy:.3g} deg , gamma = {gamma_cauchy:.3g} deg')
        print(f'Maximum pdf function M={MxsCM:.3g} mb/rad, Sampling Efficiency={eta_rejection*100:.3g}%')
    
        reaction_config = {
            'mode': reaction_data['mode'],
            'interpolator': fxsCM,
            'range_theta': [np.min(tCM), np.max(tCM)],
            'XS_max': MxsCM,
            'cauchy_param': [theta0_cauchy, gamma_cauchy],
            'masses': [Mp, Me, mt, mr],
            'atomic_masses': [A_p, A_e],
            'ions': [ion_p, ion_e]
        }
    
        # Plot distribution
        fig,ax = plt.subplots(figsize=(6, 5))
        ax.plot(tCM,xsCM,'-')
        ax.set_xlabel(r'$\theta_{CM}$ (deg)')
        ax.set_xlim([0,180])
        ax.set_xticks(np.arange(0,181,30))
        ax.set_ylabel(r'$\dfrac{d\sigma}{d\Omega}$ (mb/sr)')
        ax.set_yscale('log')
        ax.set_title(f"Total cross section: {xs_tot:.4g} mb")
        # plt.show(block=False)
        experiment_name = f'{A_p}{ion_p}_{A_e}{ion_e}'
        plt.savefig(f'./outputs/{experiment_name}_XS.png')
    
    elif reaction_data['mode'] == 'lab_data_uncorrelated':
        A_p, ion_p, Mp = ion_identifier(reaction_data['projectile'])
        if 'mass_p' in reaction_data:
            Mp = reaction_data['mass_p']
            
        A_e, ion_e, Me = ion_identifier(reaction_data['ejectile'])
        if 'mass_e' in reaction_data:
            Me = reaction_data['mass_e']
    
        Mp = Mp*amu 
        Me = Me*amu
    
        if 'data_XS' in reaction_data:
            t, xs = np.loadtxt(reaction_data['data_XS'], unpack=True)
            i_sort = np.argsort(t)
            t = t[i_sort]
            xs = xs[i_sort]
        elif 'gaussfit_XS' in reaction_data:
            A  = reaction_data.get('gaussfit_XS')['A']
            x0 = reaction_data.get('gaussfit_XS')['x0']
            a0 = reaction_data.get('gaussfit_XS')['a0']
            range_theta = reaction_data.get('gaussfit_XS')['range_theta']
            t = np.linspace(range_theta[0], range_theta[1], 901)
            xs = gaussfit(t,A,x0,a0)
        else:
            raise ValueError('Missing data on Lab Angle spectrum')
        fxs = interp1d(t, xs*np.sin(np.radians(t)), kind='linear')
    
        # Check total cross-section
        xx = np.linspace(np.min(t),np.max(t),2001) # sampling for integral
        xs_tot = integrate.trapezoid(fxs(xx), x=np.radians(xx) )*2*np.pi
        print(f"Calculate Total cross section: {xs_tot:.4g} mb")
    
        area_xs = integrate.trapezoid(fxs(xx), x=xx )
        theta0_cauchy, gammatheta_cauchy, Mxs, eta_rejection = optimize_cauchy_envelope(fxs, np.min(t), np.max(t), area_xs)
        print('Optimize rejection sampling efficiency for Angle')
        print(f'Cauchy parameters: theta_0 = {theta0_cauchy:.3g} deg , gamma = {gammatheta_cauchy:.3g} deg')
        print(f'Maximum pdf function M={Mxs:.3g}, Sampling Efficiency={eta_rejection*100:.3g}%')
    
    
        if 'data_E' in reaction_data:
            E, gE = np.loadtxt(reaction_data['data_E'], unpack=True)
            i_sort = np.argsort(E)
            E = E[i_sort]
            gE = gE[i_sort]
        elif 'gaussfit_E' in reaction_data:
            A  = reaction_data.get('gaussfit_E')['A']
            x0 = reaction_data.get('gaussfit_E')['x0']
            a0 = reaction_data.get('gaussfit_E')['a0']
            min_E, max_E = np.min(np.array(x0)-4*np.array(a0)), np.max(np.array(x0)+4*np.array(a0))
            E = np.linspace(min_E, max_E, 1001)
            gE = gaussfit(E,A,x0,a0)
        else:
            raise ValueError('Missing data on Lab Energy spectrum')
        fE = interp1d(E, gE, kind='linear')

        if 'Ep_ref' in reaction_data:
            Ep_ref = reaction_data['Ep_ref']
        else:
            raise ValueError('Missing reference energy of the ion projectile: this needs to be included to shift for other energies.')
    
        # Check energy pdf normalization
        xx2 = np.linspace(np.min(E), np.max(E), 2001)
        area_fE = integrate.trapezoid(fE(xx2), x=xx2 )    
        print(f"Energy PDF normalization: {area_fE:.5f}")
    
        E0_cauchy, gammaE_cauchy, ME, eta_rejection = optimize_cauchy_envelope(fE, np.min(E), np.max(E), area_fE)
        print('Optimize rejection sampling efficiency for Energy')
        print(f'Cauchy parameters: E_0 = {E0_cauchy:.5g} MeV , gamma = {gammaE_cauchy:.3g} MeV')
        print(f'Maximum pdf function M={ME:.3g}, Sampling Efficiency={eta_rejection*100:.3g}%')
    
        # Plot data distributions
        fig,ax = plt.subplots(2,1,figsize=(6, 5))
        ax[0].plot(xx,fxs(xx),'-')
        ax[0].set_xlabel(r'$\theta_{LAB}$ (deg)')
        ax[0].set_ylabel(r'$\dfrac{d\sigma}{d\Omega}$ (mb/sr)')
        ax[1].plot(xx2, fE(xx2),'-')
        ax[1].set_xlabel(r'$E_{LAB}$ (MeV)')
        ax[1].set_ylabel(r'p.d.f.$(E_{LAB}) (a.u.)$')
        plt.tight_layout()
        # plt.show(block=False)
        experiment_name = f'{A_p}{ion_p}_{A_e}{ion_e}'
        plt.savefig(f'./outputs/{experiment_name}_XS.png')
    
        reaction_config = {
            'mode': reaction_data['mode'],
            'interpolator': [fxs, fE],
            'range_theta': [np.min(t),np.max(t)],
            'XS_max': Mxs,
            'cauchy_param_XS': [theta0_cauchy, gammatheta_cauchy],
            'range_E': [np.min(E),np.max(E)],
            'fE_max': ME,
            'cauchy_param_E': [E0_cauchy, gammaE_cauchy],
            'Ep_ref': Ep_ref,
            'masses': [Mp, Me],
            'atomic_masses': [A_p, A_e],
            'ions': [ion_p, ion_e]
        }


    elif reaction_data['mode'] == 'lab_data':
        # 1. Identify particles and masses
        A_p, ion_p, Mp = ion_identifier(reaction_data['projectile'])
        if 'mass_p' in reaction_data:
            Mp = reaction_data['mass_p']
            
        A_e, ion_e, Me = ion_identifier(reaction_data['ejectile'])
        if 'mass_e' in reaction_data:
            Me = reaction_data['mass_e']

        Mp = Mp * amu 
        Me = Me * amu

        # 2. Reference Energy
        if 'Ep_ref' in reaction_data:
            Ep_ref = reaction_data['Ep_ref']
        else:
            raise ValueError('Missing reference energy of the ion projectile: this needs to be included to shift for other energies.')

        # 3. Build 2D Interpolator from file
        if 'data_XS_E' in reaction_data:
            # Unpack the 5 variables, including the new total cross-section
            interp_2d, range_theta, range_E, total_volume, xs_tot = create_2d_interpolator_from_file(reaction_data['data_XS_E'])
        else:
            raise ValueError('Missing 2D data file (theta, E, probability). Expected key: "data_XS_E" in reaction_data.')

        # Print the physically meaningful cross-section
        print(f"Calculated Total cross section: {xs_tot:.4g} mb")

        # 4. Optimize the 2D Cauchy envelope for Rejection Sampling
        centers, gammas, M_max, eta_rejection = optimize_cauchy_2d_envelope(
            interp_2d, range_theta, range_E, total_volume
        )
        
        theta0_cauchy, E0_cauchy = centers
        gammatheta_cauchy, gammaE_cauchy = gammas

        print('Optimize 2D rejection sampling efficiency')
        print(f'Cauchy centers: theta_0 = {theta0_cauchy:.3g} deg, E_0 = {E0_cauchy:.5g} MeV')
        print(f'Cauchy gammas:  gamma_theta = {gammatheta_cauchy:.3g} deg, gamma_E = {gammaE_cauchy:.3g} MeV')
        print(f'Maximum pdf function M = {M_max:.3g}, Sampling Efficiency = {eta_rejection*100:.3g}%')

        # 5. Plot 2D distribution (Color mesh map)
        # Create a mesh grid to evaluate the interpolator for visualization
        t_plot = np.linspace(range_theta[0], range_theta[1], 100)
        e_plot = np.linspace(range_E[0], range_E[1], 100)
        T_mesh, E_mesh = np.meshgrid(t_plot, e_plot)
        
        # Format the points as an (N, 2) array for the RegularGridInterpolator
        points_to_eval = np.column_stack((T_mesh.ravel(), E_mesh.ravel()))
        Z_plot = interp_2d(points_to_eval).reshape(T_mesh.shape)

        fig, ax = plt.subplots(figsize=(7, 5))
        c = ax.pcolormesh(T_mesh, E_mesh, Z_plot, shading='auto', cmap='viridis')
        fig.colorbar(c, ax=ax, label='Probability (a.u.)')
        ax.set_xlabel(r'$\theta_{LAB}$ (deg)')
        ax.set_ylabel(r'$E_{LAB}$ (MeV)')
        ax.set_title(f'Total cross section: {xs_tot:.4g} mb')
        plt.tight_layout()
        
        # Save the plot
        experiment_name = f'{A_p}{ion_p}_{A_e}{ion_e}'
        plt.savefig(f'./outputs/{experiment_name}_XS_E_LAB.png')
        
        # 6. Build the updated reaction configuration dictionary
        reaction_config = {
            'mode': reaction_data['mode'],
            'interpolator_2d': interp_2d,
            'range_theta': range_theta,
            'range_E': range_E,
            'M_max': M_max,
            'cauchy_param_2d': [centers, gammas],
            'Ep_ref': Ep_ref,
            'masses': [Mp, Me],
            'atomic_masses': [A_p, A_e],
            'ions': [ion_p, ion_e]
        }
    
    
    else:
        raise ValueError('Reaction type is not valid')

    return reaction_config

# ------------------------------

def build_covariance_matrix(config):
    """Builds 6x6 covariance matrix from beam config dictionary"""
    
    if config["mode"] == "covariance_6d":
        return np.array(config["matrix"])
        
    elif config["mode"] == "subspaces":
        sigma = np.zeros((6, 6))
        
        sigma[0:2, 0:2] = parse_subspace(config["x"])
        sigma[2:4, 2:4] = parse_subspace(config["y"])
        sigma[4:6, 4:6] = parse_subspace(config["z"])
        
        return sigma
        
    else:
        raise ValueError(f"Undefined beam configuration mode: {config['mode']}")

def parse_subspace(sub_config):
    """Evaluate the 2x2 sub-matrix for the specific phase subspace."""
    sub_sigma = np.zeros((2, 2))
    
    if sub_config["mode"] == "rms":
        # Uncorrelated, non-diagonal elements are zero
        sub_sigma[0, 0] = sub_config["rms_pos"]**2
        sub_sigma[1, 1] = sub_config["rms_div"]**2
        
    elif sub_config["mode"] == "twiss":
        eps = sub_config["emit"]
        beta = sub_config["beta"]
        alpha = sub_config["alpha"]
        gamma = (1 + alpha**2) / beta
        
        sub_sigma[0, 0] = eps * beta
        sub_sigma[1, 1] = eps * gamma
        sub_sigma[0, 1] = -eps * alpha
        sub_sigma[1, 0] = -eps * alpha
        
    else:
        raise ValueError(f"Undefined subspace config mode: {sub_config['mode']}")
        
    return sub_sigma

# ---------------------------------

def initialize_beam(N, beam_config, Mp, A_p, start_idx=0):
    """
    Initializes the beam either from a covariance matrix (RMS/Twiss) 
    or by loading data from an external text file.
    """
    
    # --- NEW LOGIC: Load from file with Direction Cosines ---
    if beam_config.get('mode') == 'file_data':
        if 'file_data_cache' not in beam_config:
            filename = beam_config.get('filename', 'beam.txt')
            # Columns: x[m], nx, y[m], ny, t[s], E_kin[MeV]
            beam_config['file_data_cache'] = np.loadtxt(filename)
        
        full_data = beam_config['file_data_cache']
        total_rows = len(full_data)
        
        end_idx = min(start_idx + N, total_rows)
        batch_data = full_data[start_idx:end_idx]
        actual_N = len(batch_data)
        
        if actual_N == 0:
            raise EOFError("No more particles available in the input file.")

        # 1. Extract Positions (x, y, 0)
        x0 = np.column_stack((batch_data[:, 0], batch_data[:, 2], np.zeros(actual_N)))
        
        # 2. Extract Direction Cosines (nx, ny)
        nx = batch_data[:, 1]
        ny = batch_data[:, 3]
        
        # 3. Calculate nz = sqrt(1 - nx^2 - ny^2)
        # We use np.clip to ensure the sum doesn't slightly exceed 1.0 due to precision
        sum_sq = np.clip(nx**2 + ny**2, 0, 1.0)
        nz = np.sqrt(1.0 - sum_sq)
        
        # 4. Construct p0 (Already normalized by definition)
        p0 = np.column_stack((nx, ny, nz))
        
        # 5. Time and Energy
        t0 = batch_data[:, 4]
        E0 = batch_data[:, 5]
        
        return x0, p0, t0, E0

    # --- ORIGINAL LOGIC: Covariance Matrix ---
    if 'E0/A' in beam_config:
        Ep = beam_config['E0/A'] * A_p
        pc0 = np.sqrt(Ep * (Ep + 2 * Mp))
    elif 'E0' in beam_config:
        Ep = beam_config['E0']
        pc0 = np.sqrt(Ep * (Ep + 2 * Mp))
    elif 'pc0' in beam_config:
        pc0 = beam_config['pc0']
    else:
        raise ValueError('Undefined beam energy or momentum in beam_config.')
    
    covariance_matrix = build_covariance_matrix(beam_config)
    
    # Generate random data from 6D distribution
    data = np.random.multivariate_normal(np.zeros(6), covariance_matrix, N)
    
    # Initial position
    x0 = np.column_stack((data[:, 0], data[:, 2], np.zeros(N)))
    
    # Initial direction (normalized)
    p0 = np.column_stack((data[:, 1], data[:, 3], np.ones(N)))
    p0_norm = np.linalg.norm(p0, axis=1, keepdims=True)
    p0 = p0 / p0_norm
    
    # Time distribution
    t0 = data[:, 4]
    
    # Initial momentum and energy distribution
    pc = pc0 * (1 + data[:, 5])
    E0 = np.sqrt(Mp**2 + pc**2) - Mp
    
    return x0, p0, t0, E0