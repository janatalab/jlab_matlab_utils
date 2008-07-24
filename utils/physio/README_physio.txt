README_physio.txt

The scripts and functions in this directory are an attempt to have a
unified utility for reading, parsing, and processing external
physiological data recorded in the course of fMRI or other
experiments.  Given that data recorded by different systems differs in
format, separate modules are being written to handle the data from
each system.

Important scripts
init_physmon_params.m - initializes a structure that specifies various
control parameters for parsing and processing the data. It takes an
optional input argument that populates fields that are specific to a
particular acquisition system.


Systems currently supported:
1. Dartmouth Brain Imaging Center - Data recorded using the Keithley
acquisition system.
2. UC Davis, Imaging Research Center, MATE system on Siemens 3T

