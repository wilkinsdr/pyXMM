# pyXtract
A Python interface to astronomical X-ray data reductio pipelines.

Data reduction and extraction of spectral products is handled by the various Extractor classes for different instruments and missions:

XMM-Newton:
- pyxtract.epic.EPICExtractor - for the EPIC MOS and pn CCD imaging spectrometers
- pyxtract.rgs.RGSExtractor - for the Reflection Grating Spectrometers
- pyxtract.om.OMExtractor - for the Optical Monitor

NuSTAR:
- pyxtract.nustar.NustarExtractor - for NuSTAR FPMA and FPMB

NICER:
- pyxtract.nicer.NicerExtarctor - for NICER

XRISM (coming soon):
- pyxtract.resolve.ResolveExtractor - for the XRISM Resolve microcalorimeter spectrometer
- pyxtract.xtend.XtendExtractor - for the XRISM Xtend CCDs

This automates typical data reduction tasks and the extraction of spectra and light curves and takes care of SAS environment variables, event lists, region files, etc.

pyXMM requires HEASOFT initialised in your shell environment, in addition to SAS for the XMM-Newton routines.

Initialise an EPICExtractor object for the OBSID directory and speficy either 'pn' or 'mos' as the instrument.

ex = EPICExtractor('012345678', instrument='pn', run_reduction=True)

If run_reduction is True, the standard processing will be run for a new observation (building the calibration index, ingesting the ODFs and filtering the event lists).

Then, you can start extracting data products, e.g.

ex.get_spectrum()
