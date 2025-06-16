# Ilastik Prediction Scripts

We predict object boundaries and semantic class probabilities with the [ilastik autocontext workflow](https://www.ilastik.org/documentation/autocontext/autocontext).
Autocontext performs multiple rounds of semantic pixel classification with a Random Forest Classifier based on convolutional features.
The input for each stage consists of raw data AND the predictions of the previous stage.

Unfortunately, the implementation of autocontext in ilastik does not scale to data of our size, so we emulate the workflow
and manually stack raw data and predictions.

To run all predictions stages, you will need to run the scripts in this folder in the following order:
```
./2_run_prediction.py 0
./1_make_input_data.py 1
./2_run_prediction.py 1
./1_make_input_data.py 2
./2_run_prediction.py 2
```

Note that you will need to adjust the following paths in the scripts:
- path to n5-dataset with raw data
- path to the ilastik projects
- path to python interpreter of conda environment with all dependencies
