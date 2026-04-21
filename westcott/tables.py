import numpy as np
np.set_printoptions(legacy='1.25')
from scipy.integrate import trapezoid
import pandas as pd
pd.set_option('display.max_rows', None)
import os, glob
import re
import csv

from .log_handlers import *

class CrossSectionData(object):
    __doc__="""Class to handle neutron-capture cross section data tables from 
    ENDF-B/VIII.1."""

    _ROOT = os.path.abspath(os.path.dirname(__file__))

    def __init__(self):
        from . import get_data
        self.capture_data_path = get_data('data_capture')
        capture_data_list = [x for x in glob.glob("{0}/*.csv".format(self.capture_data_path))]
        self.capture_data_dict = {}
        for c in capture_data_list:
            c_file = c.split('data_capture/')[1]
            target = c_file.split('n-capture-')[1].split('.csv')[0]
            self.capture_data_dict.update({target: c_file})
        self.capture_data_dict = dict(sorted(self.capture_data_dict.items()))

    def find_targets(self):
        """Find list of targets containing neutron-capture cross-section data.

        Arguments:
            None

        Returns:
            List object containing string representations of neutron-capture 
            target nuclides.

        Raises:
            Additional positional arguments raises a TypeError exception.

        Example:
            To find list of all available targets:
            > find_targets()
        """
        return [target for (target, value) in self.capture_data_dict.items()]

    def get_MT102(self, target):
        """Retrieve point-wise cross section data for defined target from the 
        ENDF/B-VIII [1] library: energy [eV]; cross section [b].

        Notes:
            [1] D. Brown at al., ENDF/B-VIII.0}: The 8th Major Release of the 
                Nuclear Reaction Data Library with CIELO-project Cross 
                Sections, New Standards and Thermal Scattering Data, Nucl. Data
                Sheets 148, 1 (2018).

        Arguments:
            target: A string argument that defines the target nucleus according 
                    to chemical symbol and mass.

        Returns:
            DataFrame object containing energy [eV] and cross-section [b] data 
            corresponding to target nucleus argument.

        Raises:
            Fewer or more positional arguments raises a TypeError exception.

        Example:
            To extract ENDF neutron-capture cross-section data for the 'Gd157' 
            target:
            > get_MT102('Gd157')
        """
        self.target = target
        df = None
        for k,v in self.capture_data_dict.items():
            if k==self.target:
                df = pd.read_csv("{0}/{1}".format(self.capture_data_path, v))
        if df is None:
            logger.error("No capture-gamma cross section data for target nucleus: {0}".format(target))
            return
        else:
            return df

    def sigma_ENDF(self, target):
        """Convert ENDF/B-VIII [1] energy [eV] and cross-section [b] data from 
        DataFrame object into numpy arrays allowing for interpolation and 
        integration.

        Notes:
            [1] D. Brown at al., ENDF/B-VIII.0}: The 8th Major Release of the 
                Nuclear Reaction Data Library with CIELO-project Cross 
                Sections, New Standards and Thermal Scattering Data, Nucl. Data
                Sheets 148, 1 (2018).

        Arguments:
            target: A string argument that defines the target nucleus according
                    to chemical symbol and mass.

        Returns:
            Tuple object containing two elements corresponding to target 
            nucleus argument:
        
            [0]: Numpy array of the energy in units of [eV].
            [1]: Numpy array of the cross section in units of [b].

        Raises:
            Fewer or more positional arguments raises a TypeError exception.

        Example:
            To create numpy arrays for the energy and cross-section data from 
            ENDF for the 'Gd157' target:
            > energy, cross_section = sigma_ENDF('Gd157')
        """
        self.target = target

        df = CrossSectionData.get_MT102(self,self.target)

        endf_data = df.to_numpy()
        En = endf_data[:,0]
        sigma = endf_data[:,1]

        return (En, sigma)

class NeutronFlux(CrossSectionData):
    __doc__="""Class to handle experimental neutron-flux spectra from the 
    Budapest Research Reactor (BRR) and Garching (FMR II)."""
   
    _ROOT = os.path.abspath(os.path.dirname(__file__))
   
    def __init__(self):
        from . import get_data
        self.flux_data_path = get_data('data_spectra')
        flux_data_list = [x for x in glob.glob("{0}/*.csv".format(self.flux_data_path))]
        self.flux_data_dict = {}
        for i,f in enumerate(flux_data_list):
            f_file = f.split('data_spectra/')[1]
            self.flux_data_dict.update({i: f_file})
        self.flux_data_dict = dict(sorted(self.flux_data_dict.items()))

    def find_flux(self):
        """Find available experimental neutron-flux spectra from the different
        facilities.

        Arguments:
            None

        Returns:
            List object containing two-element tuples corresponding to the 
            different experimental neutron-flux spectra:
            
            [0]: Integer (0-3) identifying experimental neutron-flux spectrum.
            [1]: String representing name of corresponding CSV-formatted file.
            
        Raises:
            Additional positional arguments raises a TypeError exception.

        Example:
            To find list of all available experimental neutron-flux spectra:
            > find_flux()
        """
        return [(i,f) for (i, f) in self.flux_data_dict.items()]

    def get_flux_df(self,flux):
        """Retrieve point-wise data for desired experimental neutron-flux 
        spectrum.

        Arguments:
            flux: An integer argument corresponding to the specified 
                  experimental neutron-flux spectrum:

                  0: Budapest Research Reactor, cold (2002).
                  1: Budapest Research Reactor, cold (2012).
                  2: Budapest Research Reactor, thermal (2002).
                  3: FRM-II Reactor, cold (2008).

        Returns:
            DataFrame object containing energy [eV] and flux [pdf] data 
            corresponding to flux argument.

        Raises:
            Fewer or more positional arguments raises a TypeError exception.

        Examples:
            To get the flux DataFrame for the cold Budapest spectrum (2002):
            > get_flux_df(0)
            To get the flux DataFrame for the cold Budapest spectrum (2012):
            > get_flux_df(1)
            To get the flux DataFrame for the thermal Budapest spectrum (2002):
            > get_flux_df(2)
            To get the flux DataFrame for the cold FRM-II spectrum (2008):
            > get_flux_df(3)
        """
        self.flux = flux
        df = None
        for k,v in self.flux_data_dict.items():
            if k==self.flux:
                df = pd.read_csv("{0}/{1}".format(self.flux_data_path, v))
        if df is None:
            logger.warning("Spectrum not defined for argument:".format(self.flux))
            return
        else:
            return df

    def get_flux(self, flux):
        """Convert experimental energy [eV] and flux [pdf] data from DataFrame 
        object into numpy arrays for allowing for interpolation and 
        integration.

        Arguments:
            flux: An integer argument corresponding to the specified 
                  experimental neutron-flux spectrum:

                  0: Budapest Research Reactor, cold (2002).
                  1: Budapest Research Reactor, cold (2012).
                  2: Budapest Research Reactor, thermal (2002).
                  3: FRM-II Reactor, cold (2008).

        Returns:
            Tuple object containing two elements corresponding to neutron-flux 
            argument:
        
            [0]: Numpy array of the energy (E) in units of [eV].
            [1]: Numpy array of the flux (dn/dE) in units of [pdf].

        Raises:
            Fewer or more positional arguments raises a TypeError exception.

        Examples:
            To get the flux data for the cold Budapest spectrum (2002):
            > energy, flux = get_flux(0)
            To get the flux data for the cold Budapest spectrum (2012):
            > energy, flux = get_flux(1)
            To get the flux data for the thermal Budapest spectrum (2002):
            > energy, flux = get_flux(2)
            To get the flux data for the cold FRM-II spectrum (2008):
            > energy, flux = get_flux(3)
        """
        self.flux = flux

        df = NeutronFlux.get_flux_df(self, self.flux)

        flux_spectrum = df.to_numpy()
        En = flux_spectrum[:,0]
        sigma = flux_spectrum[:,1]

        return (En, sigma)

class ResonanceData(NeutronFlux):
    __doc__="""Class to handle the Breit-Wigner and Reich-Moore resonances of 
    the ENDF-B/VIII.1 library."""

    _ROOT = os.path.abspath(os.path.dirname(__file__))
    
    def __init__(self):
        from . import get_data
        self.res_data_path = get_data('data_resonances')

        # Load Breit-Wigner resonances
        breit_wigner_list = [x for x in glob.glob("{0}/BreitWigner/*.csv".format(self.res_data_path))]
        breit_wigner_dict = {}
        for bw in breit_wigner_list:
            bw_file = bw.split('BreitWigner/')[1]
            target = bw_file.split('n-res-')[1].split('.csv')[0]
            breit_wigner_dict.update({target: [bw_file, 'BreitWigner']})
        #Sort dictionary by key
        bw_sorted_dict = dict(sorted(breit_wigner_dict.items()))

        # Load Reich-Moore resonances
        reich_moore_list = [x for x in glob.glob("{0}/ReichMoore/*.csv".format(self.res_data_path))]
        reich_moore_dict = {}
        for rm in reich_moore_list:
            rm_file = rm.split('ReichMoore/')[1]
            target = rm_file.split('n-res-')[1].split('.csv')[0]
            reich_moore_dict.update({target: [rm_file, 'ReichMoore']})
        #Sort dictionary by key
        rm_sorted_dict = dict(sorted(reich_moore_dict.items()))

        # Concatenate and sort the Breit-Wigner and Reich-Moore dictionaries
        merged_res_dict = {**bw_sorted_dict, **rm_sorted_dict}
        self.res_sorted_dict = dict(sorted(merged_res_dict.items()))

    def find_resonances(self,**kwargs):
        """Find list of all 'target+n' systems with resonance parameters in 
        the ENDF-B/VIII.1 [2] library: includes Breit-Wigner and Reich-Moore.

        Notes:
            [2] G. Nobre, D. Brown, R. Arcilla, R. Coles, B. Shu, Progress 
                towards the ENDF/B-VIII.1 release, Eur. Phys. J. (Web of 
                Conf.), 294, 04004 (2024).

        Arguments:
            None: All resonance parameters (Breit-Wigner and Reich-Moore).
            kwargs: Optional keyword used to specify resonance-parameter type:
                    `res='BW'` for Breit-Wigner
                    `res='RM'` for Reich-Moore

        Returns:
            List object containing target nuclides that have associated 
            resonance parameters of the form Breit-Wigner or Reich-Moore.

        Raises:
            A TypeError exception gets raised if a string argument gets passed 
            without a keyword.

        Example:
            Find all resonances:
            >find_resonances()
            Find Breit-Wigner resonances only:
            >find_resonances(res='BW')
            Find Reich-Moore resonances only:
            >find_resonances(res='RM')
        """
        self.kwargs = kwargs

        if self.kwargs == {} or self.kwargs is None:
            # return all targets
            return [target for (target, value) in self.res_sorted_dict.items()]
        else:
            for key in self.kwargs.keys():
                if key == 'res':
                    for res in kwargs.values():
                        if res.upper() == 'BW':
                            return [target for (target, value) in self.res_sorted_dict.items() if value[1]=='BreitWigner']
                        elif res.upper() == 'RM':
                            return [target for (target, value) in self.res_sorted_dict.items() if value[1]=='ReichMoore']
                        else:
                            logger.warning("Unknown keyword argument for resonance parametrizations.\nUse one of the following methods:\n`find_resonances()`\n`find_resonances(res='BW')`\n`find_resonances(res='RM')`")
                            return
                else:
                    logger.warning("Unkown key: use one of the following methods:\n`find_resonances()`\n`find_resonances(res='BW')`\n`find_resonances(res='RM')`")
                    return

    def get_res_paras(self, target):
        """Extract resonance parameters from ENDF-B/VIII.1 [2] for a defined 
        target nucleus.  The parameters are of the form Breit-Wigner or 
        Reich-Moore.

        Notes:
            [2] G. Nobre, D. Brown, R. Arcilla, R. Coles, B. Shu, Progress 
                towards the ENDF/B-VIII.1 release, Eur. Phys. J. (Web of 
                Conf.), 294, 04004 (2024).

        Arguments:
            target: A string argument that defines the target nucleus according
                    to chemical symbol and mass.
        
        Returns:
            DataFrame object containing different columns depending on the 
            resonance parameter selection for the defined target nucleus 
            argument: 

            Breit-Wigner: energy [eV]; orbital angular momentum (L); spin (J);
                          total width [eV]; neutron width [eV]; 
                          capture width [eV].

            Reich-Moore: energy [eV]; photon width [eV]; neutron width [eV]; 
                         spin (J); parity (Pi); orbital angular momentum (L).
        
        Raises:
            Fewer or more positional arguments raises a TypeError exception.

        Examples:
            Extract Breit-Wigner resonance parameters for 'Kr83' target:
            > get_res_paras('Kr83')
            Extract Reich-Moore resonance parameters for 'Gd157' target:
            > get_res_paras('Gd157')
        """
        self.target = target
        df = None
        PARAS_EXIST = False
        for nucleus, csv_object in self.res_sorted_dict.items():
            if target == nucleus and csv_object[1] == 'BreitWigner':
                PARAS_EXIST = True
                df = pd.read_csv("{0}/BreitWigner/{1}".format(self.res_data_path, csv_object[0]))
            elif target == nucleus and csv_object[1] == 'ReichMoore':
                PARAS_EXIST = True
                df = pd.read_csv("{0}/ReichMoore/{1}".format(self.res_data_path, csv_object[0]))
        if PARAS_EXIST == False:
            logger.error("No resonance parameters available for defined target or target does not exist.")
            return
                
        df_sorted = df.sort_values(by='energy')
        return df_sorted
    
