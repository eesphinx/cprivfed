# Showing the distribution bias effect

In this experiment, the alphabet size, privacy levels, and the number of datasets are fixed to 200, while the individual datasets' size varies from 50 to 400. The original data is sampled from discretized beta distributions scaled to the alphabet [0..99] and then distributed non-homogeneously between two groups of clients, each consists of 100 clients. This experiment yields the results in Figure 21. The settings can be controlled by modifying the script 'distribution_bias.sh' as follow. 

## choosing the estimation error metric
Modify the following line with the required metric
```bash
statDis="SQD"  # set to "SQD" for mean squared error, "TV" for total variation, or "EMD" for earth mover's distance
```

## no-shuffling and shuffling:
Modify the following line
```bash
shuffle_mode="shuffle" # set to "" for no shuffling,  or to "shuffle" for shuffling 
```

## non-weighted and weighted aggregation:
```bash
weights="equal" # set to "equal" for non-weighted aggregation, or "var" for weighted aggregation
```
