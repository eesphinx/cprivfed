# Varying the alphabet size

In this experiment, the privacy levels, number of datasets, and the size of each dataset are fixed, while the size of the alphabet varies from 200 to 900. The settings of the experiments can be controlled by modifying the script 'varying_alphabet.sh' as follow. 

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