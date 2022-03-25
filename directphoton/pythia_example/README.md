-------------------------------------------
Instructions for running the example shower
-------------------------------------------

1) Make sure the paths in the Makefile point to your pythia installation.

2) You need ROOT for the histogramming: 'root-config' is called for compilation.

3a) 'make example' for compiling with standard PowhegHooks

3b) 'make example_QEDQCDHooks' for compiling with separate treatment of QED and QCD radiation

4) For running, you need to specify an output root file and at least one lhe file, e.g.:
   ./example test.root pwgevents-example.lhe

Note that the pwgevents-example.lhe is for test purposes.
It contains only 1000 events, while often several millions events are needed for sufficient statistics (even with the photon radiation enhancement feature).

-------------------------------------------------
About PowhegHooks and the config file shower.conf
-------------------------------------------------
In shower.conf you can switch off features of Pythia8 such as multiparton interactions.
There are also flags concerning powheg (for a general discussion of these flags see the pythia8 online manual -> "POWHEG Merging").
What you need to know:
 - POWHEG generates additional radiation with a certain scale (stored as the SCALEUP variable in the LHE format) that would be double-counted if the parton shower would also be allowed to generate an emission above this scale. Therefore a veto for parton shower emissions above this scale is applied through Pythia's PowhegHooks.h or the modified QEDQCDPowhegHooks.h (allows the correct, separate treatment of QED and QCD scale).
 - NB: We did not find any significant difference between using the standard PowhegHooks.h or the modified QEDQCDPowhegHooks.h.
 - If you want to study the impact from vetoing, you can switch off the veto procedure with POWHEG:veto = 0.