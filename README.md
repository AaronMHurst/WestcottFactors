# Westcott $g$-factor

Python package for calculating Westcott $g$-factors, which adjust thermal-neutron capture cross sections for non-*1/v* behavior at low energies.  These are important in Neutron Activation Analysis (NAA) and Prompt Gamma-ray Activation Analysis (PGAA) when measuring nuclei that have low-energy resonances in their cross sections.  Code calculates $g$-factors for Maxwellian neutron energy distributions and individual neutron source spectra, when available.  A reference article describing this project is available on the arXiv [[1]](#1).

## Building and installation

This project can be built and installed by running the `installation.sh` script at the terminal command line of the project directory:

```Bash
$ git clone https://github.com/DMatters/WestcottFactors.git
$ cd Westcott
$ sh installation.sh
```

As an alternative the project can also be installed via `pip` since it is being concurrently maintained on the Test instance of the Python Package Index repository https://test.pypi.org/project/westcott/

```Bash
pip install -i https://test.pypi.org/simple/ westcott
```

## Testing

A suite of Python modules containing unit tests has been developed for this project.  These unit tests are located in the `tests` folder.  To run the suite and ensure they work with the local Python environment, run `tox` in the project directory where the `tox.ini` file is also located:

```Bash
$ tox -r
```

This project has the following Python-package dependencies: `numpy`, `scipy`, `pandas`, and `pytest`.  The session is automatically started after building against the required Python environment.

## Running the software

Following installation, the `westcott` scripts can be ran from any location by importing the library and making an instance of the `Westcott` class:

```Bash
$ python
```
```python
>>> import westcott
>>> gw = westcott.Westcott()
```

Various `Jupyter Notebooks` are provided in the `notebooks` folder to demonstrate workflows and methods for interacting with the functionality available to the library.

## Docstrings

Supporting docstrings for the `westcott` methods are included as part of each module.

## References
<a id="1">[1]</a>
D.A. Matters, A.M. Hurst, T. Kawano,
*"Westcott g Factors Extended to Arbitrary Neutron Energy Spectra"*,
https://doi.org/10.48550/arXiv.2602.05995
