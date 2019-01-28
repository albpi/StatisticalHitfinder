# Statistical hit-finder
# ------------------------------------------------------------------------
# Copyright 2018,  Alberto Pietrini
# Statisticalhitfinder is distributed under the terms of the Simplified BSD License.
# -------------------------------------------------------------------------

The present project refers to the journal paper by A. Pietrini et al., "A statistical approach to detect protein complexes at X-ray Free Electron Laser facilities.", Communications Physics, volume 1, Article number: 92 (2018).
The link to the paper is the following: https://www.nature.com/articles/s42005-018-0092-6.

The statistical hit-finder software is a tool specifically suited to detect protein complexes (and small biological particles in general) at XFELs (X-ray Free Electron Lasers) facilities. It is meant to be used for Flash X-ray Imaging (FXI) experiments in the peculiar case of low signal-to-noise ratio regime. It has been tested and has been shown to work at the CXI instrument at the LCLS (SLAC center, Sanford, USA). 

The script "hitfinder.py" is meant to process the raw data contained in the file cxidb-78.cxi tha can be found
in the CXIDB - Coherent Imaging Data Bank (http://www.cxidb.org/browse.html).

"hitfinder.py" executes all the steps descriped in the aforementioned paper.

In the "ipynb" folder, a Juppyter notebbok -- normalization.ipynb -- should clarify the reader how to perform the normalization of the log-likelihood scores obtained thanks the "hitfinder.py" script.

The files contained in the folder "files" have been obtained with another version of this script.
For a matter of simplicity and readability, those have been given in the cxidb-78.cxi and are used by the script "hitfinder.py" without the need of reprocessing them.
