# Varying the privacy level for half of the datasets

In this experiment, we fix the alphabet size, the number of datasets, dataset's size, and the local privacy level ε_1 for half of the datasets, while we vary the local privacy level ε_2 for the other half of the datasets from 0.1 to 3.4. The settings of the experiments can be controlled by modifying the script 'varying_leps2.sh' as follow. 

## To set the local privacy level ε_1 for the first half of datasets
```bash
leps1="0.1"   # available options are "0.1" and "3.4"
```

## To choose the estimation error metric
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