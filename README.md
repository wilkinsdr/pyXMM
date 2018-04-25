# pyXMM
A Python interface to XMM-Newton SAS.

Data reduction and extraction of spectral products is handled by the EPICExtractor class. This automates typical data reduction tasks and the extraction of spectra and light curves and takes care of SAS environment variables, event lists, region files, etc.

pyXMM requires HEASOFT and SAS initialised in your shell environment.

Initialise an EPICExtractor object for the OBSID directory and speficy either 'pn' or 'mos' as the instrument.

ex = EPICExtractor('012345678', instrument='pn', run_reduction=True)

If run_reduction is True, the standard processing will be run for a new observation (building the calibration index, ingesting the ODFs and filtering the event lists).

Then, you can start extracting data products, e.g.

ex.get_spectrum()
