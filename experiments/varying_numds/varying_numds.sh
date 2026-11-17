#!/bin/sh

make
clear
clear

# Uses two kRR mechanism: with eps1=x, eps2=y, dssize, and a list of numds values 
# 
# This script evaluates and plots the utility of central protocols LAP+EMP, GAUSS+EMP, 
# and the local protocol kRR+GIBU.  
# It uses the the executabe "federated_local.cpp" to produce data files for 
# kRR+GIBU, and "federated_central.cpp" to produce data files 
# for LAP+EMP, etc.. Then it uses the gnuplot script "varying_numds.gp" 
# to produce the distance plots. It also uses "varying_numds_time.gp" to 
# produce GIBU times plot.
#
# Utilities:
# ./federated_local <space> <sampling_mode> <statDist> <clients> <estimators> <distanceFiles>   
# ./federated_central <space> <sampling_mode> <statDist> <clients> <estimators> <datadir>   
# gnuplot -c varying_numds.gp <output_file> <stat_dist> <numds_list> <ceps1_list> <ceps2_list> <protocol_a_files label1> <protocol_b_files label2> ...
# gnuplot -c varying_numds_time.gp <output_file> <numds_list> <protocol_a_files label1> <protocol_b_files label2> ...

plot_cmd=(gnuplot -c "varying_numds.gp")
plot_time_cmd=(gnuplot -c "varying_numds_time.gp")
# Use as:  "${plot_cmd[@]}" <ARG1> <ARG2>...

# parameters
leps1=0.1
leps2=1.8 #1.8 #3.4  #10.0
delta="0.0001"
dssize=50
statDist="SQD" # EMD or SQD
sampling_mode="h" # "h" for homogeneous sampling and "nh" for non-homogeneous
numds_list="40 90 160 250 360 490 640 810" 
# numds_list="200 300 400 500 600 700 800 900" 
#numds_list="180 320 500 720 980 1280 1620 2000" 
#numds_list="200 300 500 700 1000 1300 1600 2000" 
weights="var" # "equal" or "var"
shuffle_mode="shuffle"  # Choose from:  "" and "shuffle"

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

alphabet_size=$(( hsize * vsize ))

# files to be constructed incrementally
ldp_gibu_distance_files=""
clap_n_distance_files=""
cgeom_n_distance_files=""
cgauss_n_distance_files=""
ldp_gibu_time_files=""
ceps1_list=""
ceps2_list=""

# iterate on the given list
for numds in $numds_list; do 
    # compute shuffle_numds and shuffle_dssize which are used to compute ceps1, ceps2 
    # in the shuffle mode and parameters for GIBU
    
    # Old: shuffling all datasets when they use the same eps: New: Shuffle each group separately
    # if [ $leps1 = $leps2 ]; then 
    #     shuffle_numds=1
    # else shuffle_numds=2 
    # fi
    # shuffle_numds=2

    # for kRR+GIBU
    echo " KRR eps =" ${leps1} ", KRR eps = "${leps2}
    clients="KRR ${hsize} ${vsize} ${leps1}:${dssize}:$((numds / 2)); KRR ${hsize} ${vsize} ${leps2}:${dssize}:$((numds / 2))"
    estimators="GIBU"
    ldp_gibu_distance_file=${data_dir}"/distance_GIBU_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    ldp_gibu_time_file=${data_dir}"/time_GIBU_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"

    ./federated_local "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${ldp_gibu_distance_file}" "${ldp_gibu_time_file}"

    shuffle_dssize=$(( numds * dssize / 2 ))
    if [ "$shuffle_mode" = "shuffle" ]; then
        ceps1=$(./krr_clones_feldman2020 "${leps1}" "${delta}" "${shuffle_dssize}" "${alphabet_size}")  # get amplified eps attained by ObfuscateAndCount
        ceps2=$(./krr_clones_feldman2020 "${leps2}" "${delta}" "${shuffle_dssize}" "${alphabet_size}")  # get amplified eps attained by ObfuscateAndCount
        ceps1_list=${ceps1_list}" "${ceps1}
        ceps2_list=${ceps2_list}" "${ceps2}
    else
        ceps1=${leps1}  
        ceps2=${leps2}
    fi

    # for Laplace+EMP
    echo "linear: CLAP ceps1 = "${ceps1}", ceps2= "${ceps2} 
    clients="CLAP ${hsize} ${vsize} ${ceps1}:${dssize}:$((numds / 2)); CLAP ${hsize} ${vsize} ${ceps2}:${dssize}:$((numds / 2))"
    if [ "${weights}" = "equal" ]; then
        estimators="EMPN" 
        clap_n_distance_file=${data_dir}"/distance_CLAP_EMPN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    else estimators="EMPWN"
        clap_n_distance_file=${data_dir}"/distance_CLAP_EMPWN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize"${dssize}".txt"
    fi
   ./federated_central "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${clap_n_distance_file}"
    
    # for Geom+EMP
    echo "linear: CGEOM ceps1 = "${ceps1}", ceps2= "${ceps2} 
    clients="CGEOM ${hsize} ${vsize} ${ceps1}:${dssize}:$((numds / 2)); CGEOM ${hsize} ${vsize} ${ceps2}:${dssize}:$((numds / 2))"
    if [ "${weights}" = "equal" ]; then
        estimators="EMPN" 
        cgeom_n_distance_file=${data_dir}"/distance_CGEOM_EMPN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    else estimators="EMPWN"
        cgeom_n_distance_file=${data_dir}"/distance_CGEOM_EMPWN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    fi
   ./federated_central "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${cgeom_n_distance_file}"

    # For Gauss+EMP
    echo "linear: CGAUSS ceps1 = "${ceps1}", ceps2= "${ceps2} 
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
    ldp_gibu_time_files=${ldp_gibu_time_files}" "${ldp_gibu_time_file}
done

# output
set -- ${numds_list}      # splits into positional parameters: $1=first_word $2=second_word ...
output=${space}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_from_"${1}"_to_"${8}"_dssize_"${dssize}"_"${statDist}
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
    "${plot_cmd[@]}" "${distance_output_file}" ${statDist} "${numds_list}" "${ceps1_list}" "${ceps2_list}" "${clap_n_distance_files} CLAP_EMPN" "${cgeom_n_distance_files} CGEOM_EMPN" "${cgauss_n_distance_files} CGAUSS_EMPN" "${ldp_gibu_distance_files} GIBU" 
else 
    "${plot_cmd[@]}" "${distance_output_file}" ${statDist} "${numds_list}" "${ceps1_list}" "${ceps2_list}" "${clap_n_distance_files} CLAP_EMPWN" "${cgeom_n_distance_files} CGEOM_EMPWN" "${cgauss_n_distance_files} CGAUSS_EMPWN" "${ldp_gibu_distance_files} GIBU" 
fi

# for time plot 
time_output_file=${plot_dir}"/time_"${output}
time_plot_title="alphabet size = "${alphabet_size}", dssize = "${dssize} #.", ".sprintf("{/Symbol e}_1 = %g ,{/Symbol e}_2 = %g ", real(leps1),real(leps2))
"${plot_time_cmd[@]}" "${time_output_file}" "${numds_list}" "${time_plot_title}" "${ldp_gibu_time_files} GIBU" 


echo " finished task" 



# use the following command to make the script executable ('a' for all users,'x' for executable ) 
# chmod a+x ./script.sh


# functions to return the i-th element in a range
# get_numds_at() { # called in the normal mode
#     i="$1"
#     set -- $numds_list
#     eval "echo \${$i}"
# }
#
# accessed using: $(get_numds_at $idx)
#
