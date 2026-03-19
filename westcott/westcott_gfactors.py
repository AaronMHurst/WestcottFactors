from .user import *

class Kinematics(UserSpectrum):
    __doc__="""Class to handle quantities related to neutron-beam kinematics."""

    # Constants
    eV = 1.602189e-19 #J
    kB = 1.38066e-23 #Boltzmann's constant, J/K
    m_n = 1.00866501 *1.660566e-27 #neutron mass, kg
    kB_eVK = 8.617333e-5  #Boltzmann's constant, eV/K
    h = 6.626183e-34 #Planck's constant, J*s
    h_eVs = 4.13567e-15 #Planck's constant, eV*s
   
    # Thermal quantities
    v_0 = 2200 #thermal-neutron velocity, m/s
    E_0 = 1/2 * m_n * v_0**2 / eV #eV
    T_0 = 293 # K
    lambda_0 = 1.8 # A

    def __init__(self):
        pass

    def display_thermal_properties(self) -> None:
        """Prints to screen thermal neutron properties:
        -Velocity
        -Temperature
        -Wavelength
        -Energy"""
        print(f"Velocity = {Kinematics.v_0} m/s")
        print(f"Temperature = {Kinematics.T_0} K")
        print(f"Wavelength = {Kinematics.lambda_0} A")
        print(f"Energy = {1000*Kinematics.E_0} meV")

    def display_constants(self) -> None:
        """Prints to screen various constants useful for dealing with neutron 
        beam kinematics and unit conversions."""
        print(f"1 eV = {Kinematics.eV} J")
        print(f"Neutron mass = {Kinematics.m_n} kg")
        print(f"Boltzmann constant = {Kinematics.kB} J/K")
        print(f"Boltzmann constant = {Kinematics.kB_eVK} eV/K")
        print(f"Planck constant = {Kinematics.h} J*s")
        print(f"Planck constant = {Kinematics.h_eVs} eV*s")
        
    def vel(self,En):
        """Convert neutron energy (eV) to velocity (m/s)"""
        self.En = En
        
        E_joules = self.En*Kinematics.eV
        return np.sqrt(2*E_joules/Kinematics.m_n)

    def phi_v_Maxwellian(self, T, v_array):
        """Maxwellian flux distribution at a given temperature T (K), as a function of velocity (m/s)"""
        self.T = T
        self.v_array = np.array(v_array)
        
        phi_v = []
        vt = np.sqrt(2*Kinematics.kB*self.T/Kinematics.m_n)
        for v in self.v_array:
            phi_v.append(2 * np.exp(-v**2/vt**2) * v**3/vt**4)
        return np.array(phi_v)

    def phi_v_IdealGuide(self, T, v_array):
        """Ideal Guide flux distribution at a given temperature T (K), as a function of velocity (m/s)"""
        self.T = T
        self.v_array = np.array(v_array)
        
        phi_v = []
        vt = np.sqrt(2*Kinematics.kB*self.T/Kinematics.m_n)
        for v in self.v_array:
            phi_v.append(2 * np.exp(-v**2/vt**2) * v/vt**2)
        return np.array(phi_v)

    def phi_v_arbitrary(self, En_array, phi_E_array):
        """Convert a user-defined flux distribution as a function of energy (eV) to a normalized flux distribution as a function of velocity (m/s)"""
        self.En = np.array(En_array)
        self.phi_E = np.array(phi_E_array)

        phi_v = []
        vn_array = np.sqrt(2 * self.En * Kinematics.eV / Kinematics.m_n)
        for i in range(len(vn_array)):
            phi_v.append(Kinematics.m_n * vn_array[i] * self.phi_E[i])
        return np.array(phi_v)/trapezoid(np.array(phi_v), vn_array)

    def p_v(self, v_array, phi_v_array):
        """Convert neutron flux distribution to a normalized neutron density distribution, both as a function of velocity (m/s)"""
        self.v = np.array(v_array)
        self.phi_v = np.array(phi_v_array)

        p_v = []
        for i in range(len(self.v)):
            p_v.append(self.phi_v[i] / self.v[i])
        return np.array(p_v) / trapezoid(np.array(p_v), self.v)

class Irregularity(Kinematics):
    __doc__="""Class to handle irregularity methods in the approximation of the 
    Westcott g-factor calculation."""

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    def del_0(self, v, E_resonance, Gamma):
        """Lorentzian lineshape for the irregularity function, per Molnar 
        Eqs. 1-3."""
        self.v = v
        self.E_resonance = E_resonance
        self.Gamma = Gamma
        
        K = Kinematics()
        E = 0.5 * K.m_n * self.v**2 / K.eV
        return ((self.E_resonance - K.E_0)**2 + (self.Gamma**2)/4) / ((self.E_resonance - E)**2 + (self.Gamma**2)/4)

    def gw_irregularity(self, E_resonance, Gamma, T=None, vn=np.logspace(0,5,100000)):
        """Evaluate g-factor using irregularity function method described by 
        Molnar et al. (Eqs. 1-5)."""
        self.E_resonance = E_resonance
        self.Gamma = Gamma
        self.T = T
        self.vn = vn

        K = Kinematics()
        d = []
        phi_v = K.phi_v_Maxwellian(self.T, self.vn)
        p_v = K.p_v(self.vn, phi_v)
        for i in range(len(self.vn)):
            d.append(Irregularity.del_0(self, self.vn[i], self.E_resonance, self.Gamma))
        d = np.array(d)
        return trapezoid(d * p_v, self.vn)
    
class gFactors(Irregularity):
    __doc__="""Class to handle numerical integration of complete cross-section 
    spectrum for the determination of Westcott g-factors."""

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    def gw_Maxwellian(self, T, E, sigma, vn=np.logspace(0,5,100000)):
        """Westcott g-factor according to assumed theoretical Maxwellian 
        distribution at a given neutron temperature."""
        self.T = T
        self.E = E
        self.sigma = sigma
        self.vn = vn

        K = Kinematics()
        v_sigma = K.vel(self.E)  #convert neutron energy to velocity
        phi_v = K.phi_v_Maxwellian(self.T, vn)
        p_v = K.p_v(self.vn, phi_v)
        sigma0 = np.interp(K.v_0, v_sigma, self.sigma)  #thermal cross section, barns
        sigma_interp = np.interp(self.vn, v_sigma, self.sigma)

        return 1/(sigma0 * K.v_0) * trapezoid(p_v * self.vn * sigma_interp, self.vn) / trapezoid(p_v, self.vn)

    def gw_IdealGuide(self, T, E, sigma, vn=np.logspace(0,5,100000)):
        """Westcott g-factor according to assumed theoretical Maxwellian 
        distribution at a given neutron temperature."""
        self.T = T
        self.E = E
        self.sigma = sigma
        self.vn = vn

        K = Kinematics()
        v_sigma = K.vel(self.E)  #convert neutron energy to velocity
        phi_v = K.phi_v_IdealGuide(self.T, vn)
        p_v = K.p_v(self.vn, phi_v)
        sigma0 = np.interp(K.v_0, v_sigma, self.sigma)  #thermal cross section, barns
        sigma_interp = np.interp(self.vn, v_sigma, self.sigma)

        return 1/(sigma0 * K.v_0) * trapezoid(p_v * self.vn * sigma_interp, self.vn) / trapezoid(p_v, self.vn)

    def gw_arbitrary(self, E_spectrum, phi_E_spectrum, E_endf, sigma_E_endf, vn=np.logspace(0,5,100000)):
        """Integrate to evaluate Westcott g-factor for an arbitrary neutron 
        flux distribution."""
        self.E_spectrum = E_spectrum
        self.phi_E_spectrum = phi_E_spectrum
        self.E_endf = E_endf
        self.sigma_E_endf = sigma_E_endf
        self.vn = vn

        K = Kinematics()
        sigma_0 = np.interp(K.E_0, self.E_endf, self.sigma_E_endf)
        v_spectrum = K.vel(self.E_spectrum)
        phi_v_spectrum = K.phi_v_arbitrary(self.E_spectrum, self.phi_E_spectrum)
        p_v_spectrum = K.p_v(v_spectrum, phi_v_spectrum)
        v_endf = K.vel(self.E_endf)
        sigma_v_endf = self.sigma_E_endf
        p_v_interp = np.interp(self.vn, v_spectrum, p_v_spectrum, left=0, right=0)
        sigma_interp = np.interp(self.vn, v_endf, sigma_v_endf)

        return 1/(sigma_0 * K.v_0) * trapezoid(p_v_interp * self.vn * sigma_interp, self.vn) / trapezoid(p_v_interp, self.vn)

    
