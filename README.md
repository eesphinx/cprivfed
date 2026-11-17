# Privacy-preserving Federated Distribution Estimation: Local Differential Privacy Strikes Back

This repository provides the [proofs](https://gitlab.com/esphinx/cprivFed/-/tree/main/proofs.pdf) for the theoritical results in the manuscript titled above, and a series of experiments in this manuscript that evaluate the protocols for estimating the distribution of sensitive attributes when the data points are collected by many clients untrusting each other. Every client applies a privacy mechanism to satisfy eps-DP to its local dataset. Then an aggregator collects the noisy data from all entities and estimates the global distributions. 

## Experimented methods:

1. Every local dataset is obfuscated by a central *Laplace* mechanism, then the noisy hisograms are aggregated, and normalized (LAP+N), (LAP+W+N).
2. Every local dataset is obfuscated by a central *Gauss* mechanism, then the the noisy hisograms are aggregated, and normalized (GAUSS+N), (GAUSS+W+N).
3. Every local dataset is obfuscated by a central *Geometric* mechanism, then the noisy hisograms are aggregated, and normalized (GEOM+N), (GEOM+W+N).
4. Each point in a local dataset is obfuscated by a *k-RR* mechanism, then the aggregator applies GIBU algorithm. (kRR+GIBU).

## Technical requirements

The project is implemented in C++17 and requires the following libraries, which can be installed via the [Homebrew](https://brew.sh/) utility.

- [GLPK (GNU Linear Programming Kit)](https://www.gnu.org/software/glpk/): ``` brew install glpk```

- [Armadillo (C++ library for linear algebra & scientific computing)](https://arma.sourceforge.net/): ``` brew install armadillo```

- [Gnuplot (portable command-line driven graphing utility)](http://www.gnuplot.info/): ```brew install gnuplot```

## Operating System and cross portability

The current implementation targets macOS. Porting it to another operating system requires modifying the *make* scripts in each experiment directory to ensure compatibility with the new environment.

## Experiments 

Four core experiments are located in the [experiments](https://gitlab.com/esphinx/cprivFed/-/tree/main/experiments) directory, with execution scripts available in their respective subdirectories, and details of changing their settings. 

1. *varying_leps2* : evaluating the estimation performance with varying differential privacy levels. 
2. *varying_numds* : evaluating the estimation performance with varying number of local datasets.
3. *varying_dssize* : evaluating the estimation performance with varying dataset size.
4. *varying_alphabet* : evaluating the estimation performance with varying alphabet size. 

### Adjustable Settings

The settings of each experiment can be set according to the following options.

* Estimation error metric: Mean squared error (SQD) - Total Variation (TV) - Earth Mover's Distance (EMD).
* Shuffling/no shuffling: The noisy data of the kRR may be shuffled across the clients that use the same local privacy level.
* Weighting/no weighting: The noisy histograms from clients may be aggregated with or without weighting. 
* The dataset: can be set to the Synthetic dataset, or the UCI Adult dataset.   

## Running an experiment

Consider any one of the above experiments with name <experiment_name>. The scripts of this experiment are placed in a subdirectory having the same name under [experiments](https://gitlab.com/esphinx/cprivFed/-/tree/main/experiments). 

### The estimation error
To show the estimation error in the setting of experiment <experiment_name>, execute the following steps:
```bash
cd <your_local_directory>
git clone https://gitlab.com/esphinx/cprivFed.git
cd cprivFed/experiments/<experiment_name>
./<experiment_name>.sh
```
The settings of <experiment_name> can be controlled, and the details of changing these settings are available in the page of the experiment. 

### Showing the execution time of GIBU

The scripts of experiments 2, 3, 4 produce the execution time of GIBU in their settings, in addition to the estimation error. 

## Datasets
Our experiments use two types of datasets as a source for the users' data: synthetic datasets generated at runtime, and the real [UCI Adult dataset](https://archive.ics.uci.edu/dataset/2/adult), which is in the directory [datasets](https://gitlab.com/esphinx/cprivFed/-/tree/main/datasets).

