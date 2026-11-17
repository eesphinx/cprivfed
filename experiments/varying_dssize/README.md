# Varying the dataset size

In this experiment, the alphabet size, privacy levels, and the number of datasets are fixed, while the datasets' size varies from 50 to 400. The settings of the experiments can be controlled by modifying the script 'varying_dssize.sh' as follow. 

## choosing the estimation error metric
Modify the following line with the required metric
```bash
statDis="SQD"  # set to "SQD" for mean squared error, "TV" for total variation, or "EMD" for earth mover's distance
```

## no-shuffling and shuffling:
Modify the following line
```bash
shuffle_mode="" # set to "" for no shuffling,  or to "shuffle" for shuffling 
```

## non-weighted and weighted aggregation:
```bash
weights="equal" # set to "equal" for non-weighted aggregation, or "var" for weighted aggregation
```

## To choose the dataset 
For Synthetic dataset: 

```bash
space="linear_binomial100"
hsize=1
vsize=100
```
and for the UCI Adult dataset 

```bash
space="uci_adults_ages"
hsize=1
vsize=81
```
