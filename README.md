# rsfmri-postprocessing
A MATLAB toolbox for GLM-based post-processing and analysis of resting-state fMRI data. It has been developed and extended over several years and forms the primary resting-state fMRI processing framework used across much of our group’s published work.

A key feature of the pipeline is that nuisance regression and temporal filtering can be modelled simultaneously within a single GLM, rather than applied sequentially as separate processing steps. The framework supports motion and physiological nuisance regression, CompCor, temporal filtering, spike/outlier regression, and optional generation of processed 4D fMRI data. It also includes downstream estimation of resting-state fluctuation amplitude (RSFA), ROI-to-ROI functional connectivity, and seed-to-voxel functional connectivity.

The code originated from the [GLM framework](https://github.com/MRC-CBU/riksneurotools/blob/master/GLM/glm.m) developed by Rik Henson at the MRC Cognition and Brain Sciences Unit.

See also:
[Geerligs et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5518296/).
