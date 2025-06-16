# Profiling cellular diversity in sponges informs animal cell type and nervous system evolution

## Spongilla FIBSEM Segmentation Scripts

This folder contains the scripts to segment cells, flagella and microvilli in the choanocyte chamber FIBSEM dataset.
The segmentation workflow consists of three main steps:

1. Predict object boundaries and semantic association with ilastik autocontext. The corresponding scripts can be found in `1_predict_ilastik`
2. Segment object instances with a lifted multicut solver. The corresponding scripts can be found in `2_segment_multicut`.
3. Exported the segmentation result to the [paintera](https://github.com/saalfeldlab/paintera) file format for manual correction, see `3_export`.

The folder `additional_scripts` contains the scripts to preprocess the training data for ilastik.

We have used ilastik [version 1.3.3](http://files.ilastik.org/ilastik-1.3.3b1-Linux.tar.bz2) and the `cluster_tools` software suite [version 0.1.0](https://github.com/constantinpape/cluster_tools/releases/tag/0.1.0).
To set up a python environment to use these scripts, please follow the installation instructions [here](https://github.com/constantinpape/cluster_tools#installation).
