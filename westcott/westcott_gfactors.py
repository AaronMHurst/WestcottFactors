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
        -Energy

        Arguments:
            None

        Returns:
            Void

        Raises:
            Additional positional arguments raises a TypeError exception.

        Example:
            To print thermal neutron properties:
            > display_thermal_properties()
        """
        print(f"Velocity = {Kinematics.v_0} m/s")
        print(f"Temperature = {Kinematics.T_0} K")
        print(f"Wavelength = {Kinematics.lambda_0} A")
        print(f"Energy = {1000*Kinematics.E_0} meV")

    def display_constants(self) -> None:
        """Prints to screen various constants useful for dealing with neutron 
        beam kinematics and unit conversions.

        Arguments:
            None

        Returns:
            Void

        Raises:
            Additional positional arguments raises a TypeError exception.

        Example:
            To print constants:
            > display_constants()
        """
        print(f"1 eV = {Kinematics.eV} J")
        print(f"Neutron mass = {Kinematics.m_n} kg")
        print(f"Boltzmann constant = {Kinematics.kB} J/K")
        print(f"Boltzmann constant = {Kinematics.kB_eVK} eV/K")
        print(f"Planck constant = {Kinematics.h} J*s")
        print(f"Planck constant = {Kinematics.h_eVs} eV*s")
        
    def vel(self,En):
        """Convert neutron energy to velocity: energy [eV]; velocity [m/s].

        Arguments:
            En: Neutron energy [eV] as integer or float.

        Returns:
            Floating-point value for neutron velocity [m/s].

        Raises:
            Fewer than one positional argument raises a TypeError exception.

        Example:
            To obtain neutron velocity corresponding to a neutron with an 
            energy of 5.0 eV:
            > En_5eV = vel(5.0)
        """
        self.En = En
        
        E_joules = self.En*Kinematics.eV
        return np.sqrt(2*E_joules/Kinematics.m_n)

    def phi_v_Maxwellian(self, T, v_array):
        """Evaluate Maxwellian neutron flux distribution at a given temperature,
        across an array of velocities: temperature [K]; velocity array energy
        [eV]; velocity array [m/s].

        Notes:
            [1] G. L. Molnár(Ed.),Handbook of Prompt Gamma Activation Analysis, 
            Kluwer Academic, Dordrecht, the Netherlands (2004).

        Arguments:
            T: Temperature defining Maxwellian distribution [K] as integer or
            float
            v_array: Numpy array of neutron velocities [m/s].

        Returns:
            Numpy array of Maxwellian distribution values at each velocity in
            v_array.

        Raises:
            Fewer than two positional argument raises a TypeError exception.

        Example:
            To obtain the values for a Maxwellian neutron-flux distribution, 
            defined with a temperature of 293 K, corresponding to neutron 
            velocities ranging from 1 to 10000 m/s:
            > phi_max_293K = phi_v_Maxwellian(293, np.linspace(1,10000,10000))
        """
        self.T = T
        self.v_array = np.array(v_array)
        
        phi_v = []
        vt = np.sqrt(2*Kinematics.kB*self.T/Kinematics.m_n)
        for v in self.v_array:
            phi_v.append(2 * np.exp(-v**2/vt**2) * v**3/vt**4)
        return np.array(phi_v)

    def phi_v_IdealGuide(self, T, v_array):
        """Evaluate an ideal guide neutron flux distribution at a given
        temperature, across an array of velocities: temperature [K]; 
        velocity array energy [eV]; velocity array [m/s].

        Notes:
            [1] G. L. Molnár(Ed.),Handbook of Prompt Gamma Activation Analysis,
            Kluwer Academic, Dordrecht, the Netherlands (2004).
            [2] D. A. Matters, A. M. Hurst, and T. Kawano: “Westcott g Factors
            Extended to Arbitrary Neutron Spectra,” arXiv: 2602.05995 (2026).

        Arguments:
            T: Temperature defining ideal guide distribution [K] as integer or 
            float
            v_array: Numpy array of neutron velocities [m/s].

        Returns:
            Numpy array of ideal guide distribution values at each velocity
            in v_array.

        Raises:
            Fewer than two positional argument raises a TypeError exception.

        Example:
            To obtain the values for a ideal guide neutron-flux distribution,
            defined with a temperature of 293 K, corresponding to neutron
            velocities ranging from 1 to 10000 m/s:
            > phi_ideal_293K = phi_v_IdealGuide(293, np.linspace(1,10000,10000))
        """
        self.T = T
        self.v_array = np.array(v_array)
        
        phi_v = []
        vt = np.sqrt(2*Kinematics.kB*self.T/Kinematics.m_n)
        for v in self.v_array:
            phi_v.append(2 * np.exp(-v**2/vt**2) * v/vt**2)
        return np.array(phi_v)

    def phi_v_arbitrary(self, En_array, phi_E_array):
        """Convert a user-defined neutron flux distribution as a function
        of energy to a normalized flux distribution as a function of 
        velocity: flux distribution [pdf]; energy [eV]; velocity [m/s].

        Notes:
            [1] D. A. Matters, A. M. Hurst, and T. Kawano: “Westcott g 
            Factors Extended to Arbitrary Neutron Spectra,” 
            arXiv: 2602.05995 (2026).

        Arguments:
            En_array: Numpy array of neutron energies [eV]
            phi_E_array: Numpy array of flux distribution as a function
            of energy [pdf].

        Returns:
            Numpy array of normalized flux distribution values at each
            velocity corresponding to an energy in En_array.

        Raises:
            Fewer than two positional argument raises a TypeError exception.

        Example:
            To obtain the normalized flux distribution corresponding to the
            cold-neutron source at the Budapest Research Reactor (BRR), 
            after ingesting the flux distribution from a CSV file titled
            'bnc_cold_spectrum_2012.csv' using the 'get_flux' function 
            and defining the energy array and flux array as 'En_BRR_cold'
            and 'phi_E_BRR_cold', respectively:
            > find_flux()
            > En_BRR_cold, phi_E_BRR_cold = get_flux(1)
            > phi_v_BRR_cold = phi_v_arbitrary(En_BRR_cold, phi_E_BRR_cold)
        """
        self.En = np.array(En_array)
        self.phi_E = np.array(phi_E_array)

        phi_v = []
        vn_array = np.sqrt(2 * self.En * Kinematics.eV / Kinematics.m_n)
        for i in range(len(vn_array)):
            phi_v.append(Kinematics.m_n * vn_array[i] * self.phi_E[i])
        return np.array(phi_v)/trapezoid(np.array(phi_v), vn_array)

    def p_v(self, v_array, phi_v_array):
        """Convert a neutron flux distribution as a function of velocity to
        a normalized neutron density distribution as a function of velocity:
        flux distribution [pdf]; density distribution [pdf]; velocity [m/s].
        
        Notes:
            [1] D. A. Matters, A. M. Hurst, and T. Kawano: “Westcott g
            Factors Extended to Arbitrary Neutron Spectra,” 
            arXiv: 2602.05995 (2026).

        Arguments:
            v_array: Numpy array of neutron velocities [m/s]
            phi_v_array: Numpy array of flux distribution as a function
            of velocity [pdf].

        Returns:
            Numpy array of normalized neutron density distribution values
            at each velocity in v_array.

        Raises:
            Fewer than two positional argument raises a TypeError exception.

        Example:
            To obtain the normalized density distribution corresponding to
            the cold-neutron source at the Budapest Research Reactor (BRR):,
            after ingesting the flux distribution from a CSV file titled
            'bnc_cold_spectrum_2012.csv' using the 'get_flux' function, 
            and processing it using the phi_v_arbitrary' function into a 
            flux distribution array phi_v_BRR_cold and velocity array
            v_BRR_cold:
            > find_flux()
            > En_BRR_cold, phi_E_BRR_cold = get_flux(1)
            > phi_v_BRR_cold = phi_v_arbitrary(En_BRR_cold, phi_E_BRR_cold)
            > v_BRR_cold = vel(En_BRR_cold)
            > p_v_BRR_cold = p_v(v_BRR_cold, phi_v_BRR_cold)
        """
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
        """Evaluate the neutron-capture cross section irregularity function,
        defined according to a standard Lorentzian lineshape at a given 
        resonance energy and width, at a given neutron velocity: velocity
        [m/s], resonance energy [eV], resonance width [eV].
        
        Notes:
            [1] G. L. Molnár(Ed.),Handbook of Prompt Gamma Activation 
            Analysis, Kluwer Academic, Dordrecht, the Netherlands (2004).
            (Equation 1-3).
            
        Arguments:
            v: Neutron velocity [m/s] as integer or float
            E_resonance: Energy of resonance [eV] as integer or float
            Gamma: Resonance width (FWHM) [eV] as integer or float

        Returns:
            Value of cross section irregularity function at the input
            velocity v.

        Raises:
            Fewer than three positional argument raises a TypeError exception.

        Example:
            To obtain the value of the irregularity function del_0 for
            115In (resonance energy = 1.457 eV; width = 75 meV) at a 
            neutron velocity of 2200 m/s:
            > del0_In115_thermal = del_0(2200, 1.457, 0.075)
        """
        self.v = v
        self.E_resonance = E_resonance
        self.Gamma = Gamma
        
        K = Kinematics()
        E = 0.5 * K.m_n * self.v**2 / K.eV
        return ((self.E_resonance - K.E_0)**2 + (self.Gamma**2)/4) / ((self.E_resonance - E)**2 + (self.Gamma**2)/4)

    def gw_irregularity(self, E_resonance, Gamma, T=None, vn=np.logspace(0,5,100000)):
        """Evaluate the Westcott g factor for a nucleus using the irregularity
        function method, according to a Maxwellian neutron flux distribution,
        by integrating over the neutron velocity: velocity [m/s], resonance
        energy [eV], resonance width [eV], temperature [K].
        
        Notes:
            [1] G. L. Molnár(Ed.),Handbook of Prompt Gamma Activation
            Analysis, Kluwer Academic, Dordrecht, the Netherlands (2004).
            (Equation 1-5).
            
        Arguments:
            E_resonance: Energy of cross section resonance [eV] as integer or
            float
            Gamma: Resonance width (FWHM) [eV] as integer or float
            T: Temperature defining Maxwellian distribution [K] as integer
            or float
            *vn: Numpy array of neutron velocities [m/s]. If omitted, 
            defaults to np.logspace(0,5,100000)

        Returns:
            Value of the Westcott g factor for a Maxwellian distribution 
            defined with a temperature T.

        Raises:
            Fewer than three positional argument raises a TypeError exception.

        Example:
            To obtain the value of the Westcott g factor for 115In 
            (resonance energy = 1.457 eV; width = 75 meV) at a Maxwellian
            temperature of 293 K:
            > g_In115_293K = gw_irregularity(1.457, 0.075, 293)
        """
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
        """Evaluate the Westcott g factor for a nucleus by integrating over
        the cross section in velocity space, with flux defined according 
        to a Maxwellian distribution: velocity [m/s], energy [eV], 
        cross section [b], temperature [K].
        
        Notes:
            [1] D. A. Matters, A. M. Hurst, and T. Kawano: “Westcott g
            Factors Extended to Arbitrary Neutron Spectra,” 
            arXiv: 2602.05995 (2026).
            
        Arguments:
            T: Temperature defining Maxwellian distribution [K] as integer
            or float
            E: Numpy array of neutron energies where the cross section sigma
            is defined [eV]
            sigma: Numpy array of neutron-capture cross section corresponding
            to the energies in the array E [b]
            *vn: Numpy array of neutron velocities [m/s]. If omitted, 
            defaults to np.logspace(0,5,100000)

        Returns:
            Value of the Westcott g factor for a Maxwellian distribution
            defined with a temperature T.

        Raises:
            Fewer than three positional argument raises a TypeError exception.

        Example:
            To obtain the value of the Westcott g factor for 115In, 
            using the neutron-capture cross section from ENDF/B-VIII.1
            obtained by calling the 'sigma_ENDF' function to produce 
            arrays E_sigma_In115 and sigma_In115, at a Maxwellian
            temperature of 293 K:
            > E_sigma_In115, sigma_In115 = gw.sigma_ENDF('In115')
            > g_max_293K = gw_Maxwellian(293, E_sigma_In115, sigma_In115)
        """
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
        """Evaluate the Westcott g factor for a nucleus by integrating over
        the cross section in velocity space, with flux defined according 
        to an ideal neutron guide distribution: velocity [m/s], energy [eV],
        cross section [b], temperature [K].
        
        Notes:
            [1] D. A. Matters, A. M. Hurst, and T. Kawano: “Westcott g
            Factors Extended to Arbitrary Neutron Spectra,” 
            arXiv: 2602.05995 (2026).
            
        Arguments:
            T: Temperature defining Maxwellian distribution [K] as integer
            or float
            E: Numpy array of neutron energies where the cross section 
            sigma is defined [eV]
            sigma: Numpy array of neutron-capture cross section 
            corresponding to the energies in the array E [b]
            *vn: Numpy array of neutron velocities [m/s]. If omitted,
            defaults to np.logspace(0,5,100000)

        Returns:
            Value of the Westcott g factor for an ideal guide distribution
            defined with a temperature T.

        Raises:
            Fewer than three positional argument raises a TypeError exception.

        Example:
            To obtain the value of the Westcott g factor for 115In, using
            the neutron-capture cross section from ENDF/B-VIII.1
            obtained by calling the 'sigma_ENDF' function to produce 
            arrays E_sigma_In115 and sigma_In115, at an ideal guide
            temperature of 293 K:
            > E_sigma_In115, sigma_In115 = gw.sigma_ENDF('In115')
            > g_ideal_293K = gw_IdealGuide(293, E_sigma_In115, sigma_In115)
        """
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
        """Evaluate the Westcott g factor for a nucleus by integrating over
        the cross section in velocity space, with flux defined according to
        a user-specified distribution: velocity [m/s], energy [eV], 
        cross section [b].
        
        Notes:
            [1] D. A. Matters, A. M. Hurst, and T. Kawano: “Westcott g
            Factors Extended to Arbitrary Neutron Spectra,” 
            arXiv: 2602.05995 (2026).
            
        Arguments:
            E_spectrum: Numpy array of neutron energies where the flux
            spectrum is defined [eV]
            phi_E_spectrum: Numpy array of flux distribution as a function
            of energy [pdf]
            E_endf: Numpy array of neutron energies where the cross section
            sigma_E_endf is defined [eV]
            sigma_E_endf: Numpy array of neutron-capture cross section
            corresponding to the energies in the array E_endf [b]
            *vn: Numpy array of neutron velocities [m/s]. If omitted,
            defaults to np.logspace(0,5,100000)

        Returns:
            Value of the Westcott g factor for a user-specified flux
            distribution.

        Raises:
            Fewer than four positional argument raises a TypeError exception.

        Example:
            To obtain the value of the Westcott g factor for 115In, using
            the neutron-capture cross section from ENDF/B-VIII.1
            obtained by calling the 'sigma_ENDF' function to produce arrays
            E_sigma_In115 and sigma_In115, and a flux distribution
            corresponding to the cold-neutron source at the Budapest 
            Research Reactor (BRR), after ingesting the distribution 
            from a CSV file titled 'bnc_cold_spectrum_2012.csv' using 
            the 'get_flux' function and defining the energy array and 
            flux array as 'En_BRR_cold' and 'phi_E_BRR_cold', respectively:
            > find_flux()
            > En_BRR_cold, phi_E_BRR_cold = get_flux(1)
            > E_sigma_In115, sigma_In115 = gw.sigma_ENDF('In115')
            > g_BRR_cold = gw_arbitrary(En_BRR_cold, phi_E_BRR_cold, 
                           E_sigma_In115, sigma_In115)
        """
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

    
