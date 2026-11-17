#!/bin/sh

make
clear
clear

# Uses two kRR mechanism: with leps1=x, and eps2 taking values from y to z
# 
# This script evaluates and plot the utility of central protocols LAP+EMP, GAUSS+EMP, 
# and the local protocol kRR+GIBU with varing privacy level (leps2) required by 
# the second group of clients.  
# This script uses the program "federated_local.cpp" to produce data files for 
# kRR+GIBU, and the program "federated_central.cpp" to produce data files 
# for LAP+EMP. Then it uses the gnuplot script "varying_leps2.gp" to produce
# the plots.

# Utilities:
# ./federated_local <space> <sampling_mode> <statDist> <clients> <estimators> <distanceFiles>   
# ./federated_central <space> <sampling_mode> <statDist> <clients> <estimators> <datadir>   
# gnuplot -c varying_leps2.gp "<output_file> <stat_dist> <leps2_list> <ceps2_list> "<protocol_a_files label1>" "<protocol_b_files label2>" ...

plot_cmd=(gnuplot -c "varying_leps2.gp")
# Use as:  "${plot_cmd[@]}" <ARG1> <ARG2>...

# parameters
dssize=50
numds=200
delta=0.0001
leps1=3.4 #3.4 #0.1 # 10.0
leps2_list="0.1 0.5 1.0 1.5 2.0 2.5 3.0 3.4"
# leps2_list="1.0 1.9 3.6 4.9 6.2 7.5 8.8 10.1"
# leps2_list="0.1 1.5 2.9 4.3 5.7 7.1 8.5 10.0"
statDist="SQD" # EMD or SQD
sampling_mode="h" # "h" for homogeneous sampling and "nh" for non-homogeneous
weights="var" # "equal" or "var"
shuffle_mode="shuffle"  # Choose "" or "shuffle"

# data and plot directories
data_dir="results_"${weights}"W"
if [ "${shuffle_mode}" = "shuffle" ]; then
    data_dir="${data_dir}_shuffle"
fi
plot_dir="."

# space
# space="linear_binomial100"
# hsize=1
# vsize=100
space="uci_adults_ages"
hsize=1
vsize=81
# space="sf_12km_8km"
# hsize=24
# vsize=16

alphabet_size=$((hsize * vsize))


# files to be constructed incrementally
ldp_gibu_distance_files=""
clap_n_distance_files=""
cgeom_n_distance_files=""
cgauss_n_distance_files=""
ceps2_list=""

# iterate on the given list
for leps2 in ${leps2_list}; do
    # compute shuffle_numds and shuffle_dssize which are used to computed ceps1, ceps2 
    # in the shuffle mode and parameters for GIBU
    
    # Old: shuffling all datasets when they use the same eps: New: Shuffle each group separately
    # if [ $leps1 = $leps2 ]; then 
    #     shuffle_numds=1
    # else shuffle_numds=2 
    # fi

    # for kRR+GIBU
    echo "KRR eps = " ${leps1} ", KRR eps = "${leps2}
    clients="KRR ${hsize} ${vsize} ${leps1}:${dssize}:$((numds / 2)); KRR ${hsize} ${vsize} ${leps2}:${dssize}:$((numds / 2))"
    estimators="GIBU"
    ldp_gibu_distance_file=${data_dir}"/distance_GIBU_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    
    ./federated_local "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${ldp_gibu_distance_file}"
    
    shuffle_dssize=$(( numds * dssize / 2 ))
    if [ "$shuffle_mode" = "shuffle" ]; then
        # all noisy data points are shuffled using a trusted infrastructure
        ceps1=$(./krr_clones_feldman2020 "${leps1}" "${delta}" "${shuffle_dssize}" "${alphabet_size}")  # get amplified eps1 
        ceps2=$(./krr_clones_feldman2020 "${leps2}" "${delta}" "${shuffle_dssize}" "${alphabet_size}")  # get amplified eps2
        ceps2_list=${ceps2_list}" "${ceps2}
    else
        ceps1=${leps1}
        ceps2=${leps2}
    fi

 
    # for Laplace+EMP
    clients="CLAP ${hsize} ${vsize} ${ceps1}:${dssize}:$((numds / 2)); CLAP ${hsize} ${vsize} ${ceps2}:${dssize}:$((numds / 2))"
    if [ "${weights}" = "equal" ]; then
        estimators="EMPN" 
        clap_n_distance_file=${data_dir}"/distance_CLAP_EMPN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    else estimators="EMPWN"
        clap_n_distance_file=${data_dir}"/distance_CLAP_EMPWN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    fi
   ./federated_central "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${clap_n_distance_file}"

    # for Geom+EMP
    clients="CGEOM ${hsize} ${vsize} ${ceps1}:${dssize}:$((numds / 2)); CGEOM ${hsize} ${vsize} ${ceps2}:${dssize}:$((numds / 2))"
    if [ "${weights}" = "equal" ]; then
        estimators="EMPN" 
        cgeom_n_distance_file=${data_dir}"/distance_CGEOM_EMPN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    else estimators="EMPWN"
        cgeom_n_distance_file=${data_dir}"/distance_CGEOM_EMPWN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    fi
   ./federated_central "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${cgeom_n_distance_file}"

    # For Gauss+EMP
    clients="CGAUSS ${hsize} ${vsize} ${ceps1} ${delta}:${dssize}:$((numds / 2)); CGAUSS ${hsize} ${vsize} ${ceps2} ${delta}:${dssize}:$((numds / 2))"
    if [ "${weights}" = "equal" ]; then
        estimators="EMPN" 
        cgauss_n_distance_file=${data_dir}"/distance_CGAUSS_EMPN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    else estimators="EMPWN"
        cgauss_n_distance_file=${data_dir}"/distance_CGAUSS_EMPWN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    fi
   ./federated_central "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${cgauss_n_distance_file}"

    ldp_gibu_distance_files=${ldp_gibu_distance_files}" "${ldp_gibu_distance_file}
    clap_n_distance_files=${clap_n_distance_files}" "${clap_n_distance_file}
    cgeom_n_distance_files=${cgeom_n_distance_files}" "${cgeom_n_distance_file}
    cgauss_n_distance_files=${cgauss_n_distance_files}" "${cgauss_n_distance_file}

done

# output
set -- ${leps2_list}      # splits into positional parameters: $1=first_word $2=second_word ...
output=${space}"_leps1_"${leps1}"_leps2_from_"${1}"_to_"${8}"_numds_"${numds}"_dssize_"${dssize}"_"${statDist}
if [ "${weights}" = "equal" ]; then
    output=${output}"_equalW"
else
    output=${output}"_varW"
fi
if [ "${shuffle_mode}" = "shuffle" ]; then
    output=${output}"_shuffle"
fi

# to plot distance
distance_output_file=${plot_dir}"/"${output}
if [ "${weights}" = "equal" ]; then
    "${plot_cmd[@]}" "${distance_output_file}" ${statDist} "${leps2_list}" "${ceps2_list}" "${clap_n_distance_files} CLAP_EMPN" "${cgeom_n_distance_files} CGEOM_EMPN" "${cgauss_n_distance_files} CGAUSS_EMPN" "${ldp_gibu_distance_files} GIBU" 
else 
    "${plot_cmd[@]}" "${distance_output_file}" ${statDist} "${leps2_list}" "${ceps2_list}" "${clap_n_distance_files} CLAP_EMPWN" "${cgeom_n_distance_files} CGEOM_EMPWN" "${cgauss_n_distance_files} CGAUSS_EMPWN" "${ldp_gibu_distance_files} GIBU" 
fi


echo " finished task" 



# use the following command to make the script executable ('a' for all users,'x' for executable ) 
# chmod a+x ./script.sh
