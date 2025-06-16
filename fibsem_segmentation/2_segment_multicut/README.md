# Multicut Segmentation Scripts

We segment cells, flagella and microvilli into instances by solving a Lifted Multicut problem.
Conceptually, this workflow proceeds in 3 steps:

1. Compute a watershed over-segmentation based on a height map derived from the ilastik boundary predictions.
2. Compute connected components of semantic predictions for microvilli and flagella.
3. Set up Lifted Multicut Problem:
  - Derive unary potential for edges in the Region Adjacency Graph from boundary predictions.
  - Introduce lifted edges between fragments in the over-segmentation that were mapped to the SAME 
    semantic components computed in 2. and add an attractive potential.
4. Solve Lifted Multicut problem with hierarchical solvert

Run the scripts in the following order:
```
./1_make_bmap.py              # extracts boundary map from ilastik predictions
./2_run_semantic.py           # compute connected components of the ilastik predictions for microvilli and flagella
./3_join_class_predictions.py # join connected components of microvilli and flagella into a single volume
./4_run_lmc.py                # compute over-segmentatio, lifted multicut problem and solve lifted multicut hierarchically
```

Note that you will need to adjust the following paths in the scripts:
- path to n5-dataset with raw data and ilastik predictions
- path to python interpreter of conda environment with all dependencies
