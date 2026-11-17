#!/bin/sh

make
clear
clear

# Uses two kRR mechanism: with leps1=x, leps2=y, numds, and a list of dssize values 
# 
# This script evaluates and plots the utility of central protocols LAP+EMP, GAUSS+EMP, 
# and the local protocol kRR+GIBU.  
# It uses the the executabe "federated_local.cpp" to produce data files for 
# kRR+GIBU, and "federated_central.cpp" to produce data files 
# for LAP+EMP, etc.. Then it uses the gnuplot script "varying_dssize.gp" 
# to produce the distance plots. It also uses "varying_dssize_time.gp" to 
# produce GIBU times plot.
#
# Utilities:
# ./federated_local <space> <sampling_mode> <statDist> <clients> <estimators> <distanceFiles>   
# ./federated_central <space> <sampling_mode> <statDist> <clients> <estimators> <distanceFiles>   
# gnuplot -c varying_dssize.gp <output_file> <stat_dist> <dssize_list> <ceps1_list> <ceps2_list> <protocol_a_files label1> <protocol_b_files label2> ...
# gnuplot -c varying_dssize_time.gp <output_file> <dssize_list> <protocol_a_files label1> <protocol_b_files label2> ...

plot_cmd=(gnuplot -c "varying_dssize.gp")
plot_time_cmd=(gnuplot -c "varying_dssize_time.gp")
# Use as:  "${plot_cmd[@]}" <ARG1> <ARG2>...

# parameters
leps1=0.1
leps2=3.4 #3.4  #10.0
delta="0.0001"
numds=200  
dssize_list="50 100 150 200 250 300 350 400" 
statDist="SQD" # EMD or SQD
sampling_mode="h" # "h" for homogeneous sampling and "nh2" for non-homogeneous
weights="var" # "equal" or "var"
shuffle_mode="shuffle"  # Choose from:  "" and "shuffle"

# for non-homogeneous: just set sampling_mode ="nh2" and shuffle_mode="shuffle"

# data and plot directories
data_dir="results_"${weights}"W"
if [ "${shuffle_mode}" = "shuffle" ]; then
    data_dir="${data_dir}_shuffle"
fi
plot_dir="."


# space
# space="linear_InvDiscGaussA100" # nh case 1
# space="linear_DiscBeta100" # nh case 2
# space="linear_InvDiscGaussB100" # nh case 3
# space="linear_Bimodal100" # nh case 4
space="linear_binomial100" 
hsize=1
vsize=100
# space="uci_adults_ages"
# hsize=1
# vsize=81
# space="sf_12km_8km"
# hsize=24
# vsize=16

alphabet_size=$((hsize * vsize))

# files to be constructed incrementally
ldp_gibu_distance_files=""
clap_n_distance_files=""
cgeom_n_distance_files=""
cgauss_n_distance_files=""
ldp_gibu_time_files=""
ceps1_list=""
ceps2_list=""

# iterate on the given list
for dssize in $dssize_list; do
    # for kRR + GIBU
    echo " KRR eps =" ${leps1} ", KRR eps = "${leps2}
    clients="KRR ${hsize} ${vsize} ${leps1}:${dssize}:$((numds / 2)); KRR ${hsize} ${vsize} ${leps2}:${dssize}:$((numds / 2))"
    estimators="GIBU"
    ldp_gibu_distance_file=${data_dir}"/distance_GIBU_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
    ldp_gibu_time_file=${data_dir}"/time_GIBU_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"

    ./federated_local "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${ldp_gibu_distance_file}" "${ldp_gibu_time_file}"
    
    shuffle_dssize=$(( numds * dssize / 2 )) # number of data points used to compute the amplified privacy due to shuffling
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
        clap_n_distance_file=${data_dir}"/distance_CLAP_EMPWN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_"${dssize}".txt"
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
set -- ${dssize_list}      # splits into positional parameters: $1=first_word $2=second_word ...
output=${space}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize_from_"${1}"_to_"${8}"_"${statDist}
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
    "${plot_cmd[@]}" "${distance_output_file}" ${statDist} "${dssize_list}" "${ceps1_list}" "${ceps2_list}" "${clap_n_distance_files} CLAP_EMPN" "${cgeom_n_distance_files} CGEOM_EMPN" "${cgauss_n_distance_files} CGAUSS_EMPN" "${ldp_gibu_distance_files} GIBU" 
else 
    "${plot_cmd[@]}" "${distance_output_file}" ${statDist} "${dssize_list}" "${ceps1_list}" "${ceps2_list}" "${clap_n_distance_files} CLAP_EMPWN" "${cgeom_n_distance_files} CGEOM_EMPWN" "${cgauss_n_distance_files} CGAUSS_EMPWN" "${ldp_gibu_distance_files} GIBU" 
fi

# for time plot 
time_output_file=${plot_dir}"/time_"${output}
time_plot_title="alphabet size = "${alphabet_size}", number of clients = "${numds} #.", ".sprintf("{/Symbol e}_1 = %g ,{/Symbol e}_2 = %g ", real(leps1),real(leps2))
"${plot_time_cmd[@]}" "${time_output_file}" "${dssize_list}" "${time_plot_title}" "${ldp_gibu_time_files} GIBU" 

echo " finished task" 

# use the following command to make the script executable ('a' for all users,'x' for executable ) 
# chmod a+x ./script.sh
