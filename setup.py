import setuptools

requirements = ["numpy", "scipy", "pandas"]

setuptools.setup(
    name="westcott",
    version="1.0.0",
    url="https://github.com/DMatters/WestcottFactors",
    author="David Matters and Aaron M. Hurst",
    author_email="amhurst@berkeley.edu",
    description="Calculations of Westcott g-factors.",
    long_description=open('README.md').read(),
    license_files=('LICENSE'),
    #packages=setuptools.find_packages(),
    packages=setuptools.find_namespace_packages(),
    #install_requires=[],
    install_requires=requirements,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.13',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    include_package_data=True,
    package_data={'': ['data_capture/*.csv', 'data_spectra/*.csv',
                       'data_resonances/BreitWigner/*.csv',
                       'data_resonances/ReichMoore/*.csv']},
)
