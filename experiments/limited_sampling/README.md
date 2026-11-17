# Showing the limited sampling effect

We consider 200 clients among which, the first client requires weak privacy eps1=50, and the other clients require strong (high) privacy eps2=3.0. The first (odd) client has 10 data points, while each other one has dssize2 which varies from 90 to 1000. The objective is to show the limitation of the weighting method in this case. The settings of the experiments can be controlled by modifying the script 'limited_sampling.sh' as follow. 

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
