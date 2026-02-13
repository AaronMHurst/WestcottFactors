# Westcott_g-factor

Python package for calculating Westcott $g$-factors, which adjust thermal-neutron capture cross sections for non-*1/v* behavior at low energies.  These are important in Neutron Activation Analysis (NAA) and Prompt Gamma-ray Activation Analysis (PGAA) when measuring nuclei that have low-energy resonances in their cross sections.  Code calculates $g$-factors for Maxwellian neutron energy distributions and individual neutron source spectra, when available.

## Building and installation

This project can be built and installed by running the `installation.sh` script at the terminal command line of the project directory:

```bash
$ git clone https://github.com/DMatters/WestcottFactors.git
$ cd Westcott
$ sh installation.sh
```

As an alternative the project can also be installed via `pip` since it is being concurrently maintained on the Test instance of the Python Package Index repository https://test.pypi.org/project/westcott/

```bash
pip install -i https://test.pypi.org/simple/ westcott
```
