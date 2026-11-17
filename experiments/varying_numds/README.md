# Varying the number of clients

In this experiment, the alphabet size, privacy levels, and the datasets' size are fixed, while the number of clients varies from 40 to 810. This experiment yields the results in Figures 6,7 (without weighting), Figures 13,14 (with weighting), and Figure 18(b) (GIBU execution time). The settings are controlled by modifying the script 'varying_numds.sh' as follow. 

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
