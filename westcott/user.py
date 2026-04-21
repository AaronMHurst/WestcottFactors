from .tables import *

class UserSpectrum(ResonanceData):
    __doc__="""Class to handle an arbitrary two-column CSV format neutron 
    flux spectrum."""

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    def import_spectrum(self,csv_filename):
        """Import arbitrary neutron flux spectrum as a function of energy in 
        CSV format: column 1 - energy; column 2 - flux.  It is assumed that the 
        first line of the CSV file is a header line describing the 
        comma-separated variables in each column.

        Arguments:
            csv_filename: String representation of the filename corresponding 
                          to the two-column CSV-formatted file:
 
                          Column 1: Energy
                          Column 2: Flux

        Returns:
            Tuple object containing two elements:
 
            [0]: Numpy array of pointwise energy data.
            [1]: Numpy array of pointwise flux data.

        Raises:
            Fewer or more positional arguments raises a TypeError exception.

        Example:
            For a CSV file "my_flux_data.csv" containing approproate 
            comma-separated `energy`, `flux` data, to load the contents into 
            numpy arrays:
            > energy, flux = import_spectrum("my_flux_data.csv")
        """
        self.csv_filename = csv_filename
        with open(self.csv_filename, mode='r', encoding='utf-8-sig') as file:
            csvFile = csv.reader(file)
            next(csvFile, None)  # Skip header line
            En = []
            dndE = []
            BAD_CSV_FORMAT = False
            for lines in csvFile:
                if (len(lines) < 2) or (len(lines) > 2):
                    BAD_CSV_FORMAT = True
                    logger.warning("Imported CSV file must have two columns only; tuple object will not be returned.")
                    break
                else:
                    En.append(float(lines[0]))
                    dndE.append(float(lines[1]))
            if BAD_CSV_FORMAT == False:
                return(np.array(En), np.array(dndE))

