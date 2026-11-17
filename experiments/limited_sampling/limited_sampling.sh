#!/bin/sh

make
clear
clear

# limited-sampling setting with client1 has <dssize1> points, and all others 
# have <dssize2>.
#
# Uses two kRR mechanism: with eps1=x, eps2=y, numds, and a range of dssize 
# 
# We consider 200 clients among which, the first client requires weak privacy eps1, 
# and the other clients require strong (high) privacy eps2. The first (odd) client 
# has <dssize1> data points, while each other one has <dssize2>. 
# We set <dssize1> to 10 points, and vary <dssize2> from 90 to 1000.
# The objective is to show the limitation of the weighting method in this case.
# 
# Under the above setting, this script evaluates and plot the utility 
# of the central protocols LAP+EMP, etc., and the local kRR+GIBU.  
# This script uses the the program "federated_local.cpp" to produce data files 
# for kRR+GIBU, and the program "federated_central.cpp" to produce data files 
# for LAP+EMP, GEOM+EMP, and GAUSS+EMP. 
# Then plots the data using the gnuplot script "varying_dssize.gp".

# Utilities:
# ./federated_local <space> <sampling_mode> <statDist> <clients> <estimators> <distanceFiles>   
# ./federated_central <space> <sampling_mode> <statDist> <clients> <estimators> <distanceFiles>   
# gnuplot -c limited_sampling.gp <output_file> <stat_dist> <dssize_list> <ceps1_list> <ceps2_list> <protocol_a_files label1> <protocol_b_files label2> ...

plot_cmd=(gnuplot -c "limited_sampling.gp")
# Use as:  "${plot_cmd[@]}" <ARG1> <ARG2>...

# parameters:
leps1=50.0
leps2=3.0 #3.4  #10.0
delta="0.0001"
numds=200  # this should be large
# dssize_list="5 10 15 20 25 30 35 40" 
dssize1="10"
dssize2_list="90 160 250 360 490 640 810 1000" 
# dssize2_list="50 100 150 200 250 300 350 400" 
statDist="SQD" # EMD or SQD
sampling_mode="h" # "h" for homogeneous sampling and "nh" for non-homogeneous
weights="var" # "equal" or "var"
shuffle_mode="shuffle"  # Choose from:  "" and "shuffle"


# for "nh" case: just set samling_mode ="nh" and shuffle_mode="shuffle"

# data and plot directories
data_dir="results_"${weights}"W"
if [ "$shuffle_mode" = "shuffle" ]; then
    data_dir="${data_dir}""_shuffle"
fi
plot_dir="."

# space
# space="linear_InvDiscGaussA100" # for nh case 1
# space="linear_DiscBeta100" # for nh case 2
# space="linear_InvDiscGaussB100" # for nh case 3
# space="linear_Bimodal100" # for nh case 4
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

# ceps1_list=""
ceps2_list=""

# iterate on the given list
for dssize2 in $dssize2_list; do
    # for kRR + GIBU
    echo " KRR eps = "${leps1}", KRR eps = "${leps2}", dssize2 = "${dssize2}
    freq1=1
    freq2=$(( numds - 1 ))
    mechanism1="KRR ${hsize} ${vsize} ${leps1}"
    mechanism2="KRR ${hsize} ${vsize} ${leps2}"
    clients="${mechanism1}:${dssize1}:${freq1};${mechanism2}:${dssize2}:${freq2}"
    estimators="GIBU"
    ldp_gibu_distance_file=${data_dir}"/distance_GIBU_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize1_"${dssize1}"_dssize2_"${dssize2}".txt"
    
    ./federated_local "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${ldp_gibu_distance_file}"
    
    shuffle_dssize1=$(( freq1 * dssize1 )) # number of data points used to compute the amplified privacy due to shuffling
    shuffle_dssize2=$(( freq2 * dssize2 )) # number of data points used to compute the amplified privacy due to shuffling
    if [ "$shuffle_mode" = "shuffle" ]; then
        ceps1=$(./krr_clones_feldman2020 "${leps1}" "${delta}" "${shuffle_dssize1}" "${alphabet_size}")  # get amplified eps attained by ObfuscateAndCount
        ceps2=$(./krr_clones_feldman2020 "${leps2}" "${delta}" "${shuffle_dssize2}" "${alphabet_size}")  # get amplified eps attained by ObfuscateAndCount
        # ceps1_list=${ceps1_list}" "${ceps1}
        ceps2_list=${ceps2_list}" "${ceps2}
    else
        ceps1=${leps1}
        ceps2=${leps2}
    fi

    # for CLAP+EMP
    echo "linear: CLAP ceps1 = "${ceps1}", ceps2= "${ceps2} 
    mechanism1="CLAP ${hsize} ${vsize} ${ceps1}"
    mechanism2="CLAP ${hsize} ${vsize} ${ceps2}"
    clients="${mechanism1}:${dssize1}:${freq1};${mechanism2}:${dssize2}:${freq2}"
    if [ "${weights}" = "equal" ]; then
        estimators="EMPN" 
        clap_n_distance_file=${data_dir}"/distance_CLAP_EMPN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize1_"${dssize1}"_dssize2_"${dssize2}".txt"
    else estimators="EMPWN"
        clap_n_distance_file=${data_dir}"/distance_CLAP_EMPWN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize1_"${dssize1}"_dssize2_"${dssize2}".txt"
    fi
    ./federated_central "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${clap_n_distance_file}"
    

    # for CGeom+EMP
    echo "linear: CGEOM ceps1 = "${ceps1}", ceps2= "${ceps2} 
    mechanism1="CGEOM ${hsize} ${vsize} ${ceps1}"
    mechanism2="CGEOM ${hsize} ${vsize} ${ceps2}"
    clients="${mechanism1}:${dssize1}:${freq1};${mechanism2}:${dssize2}:${freq2}"
    if [ "${weights}" = "equal" ]; then
        estimators="EMPN" 
        cgeom_n_distance_file=${data_dir}"/distance_CGEOM_EMPN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize1_"${dssize1}"_dssize2_"${dssize2}".txt"
    else estimators="EMPWN"
        cgeom_n_distance_file=${data_dir}"/distance_CGEOM_EMPWN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize1_"${dssize1}"_dssize2_"${dssize2}".txt"
    fi
   ./federated_central "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${cgeom_n_distance_file}"

    
    # For Gauss+EMP 
    echo "linear: CGAUSS ceps1 = "${ceps1}", ceps2= "${ceps2} 
    mechanism1="CGAUSS ${hsize} ${vsize} ${ceps1} ${delta}"
    mechanism2="CGAUSS ${hsize} ${vsize} ${ceps2} ${delta}"
    clients="${mechanism1}:${dssize1}:${freq1};${mechanism2}:${dssize2}:${freq2}"
    if [ "${weights}" = "equal" ]; then
        estimators="EMPN" 
        cgauss_n_distance_file=${data_dir}"/distance_CGAUSS_EMPN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize1_"${dssize1}"_dssize2_"${dssize2}".txt"
    else estimators="EMPWN"
        cgauss_n_distance_file=${data_dir}"/distance_CGAUSS_EMPWN_"${space}"_"${statDist}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize1_"${dssize1}"_dssize2_"${dssize2}".txt"
    fi
    ./federated_central "${space}" "${sampling_mode}" "${statDist}" "${clients}" "${estimators}" "${cgauss_n_distance_file}"
    

    ldp_gibu_distance_files=${ldp_gibu_distance_files}" "${ldp_gibu_distance_file}
    clap_n_distance_files=${clap_n_distance_files}" "${clap_n_distance_file}
    cgeom_n_distance_files=${cgeom_n_distance_files}" "${cgeom_n_distance_file}
    cgauss_n_distance_files=${cgauss_n_distance_files}" "${cgauss_n_distance_file}
done


######### plotting

# output
set -- ${dssize2_list}      # splits into positional parameters: $1=first_word $2=second_word ...
output=${space}"_leps1_"${leps1}"_leps2_"${leps2}"_numds_"${numds}"_dssize1_"${dssize1}"_dssize2_from_"${1}"_to_"${8}"_"${statDist}
if [ "${weights}" = "equal" ]; then
    output=${output}"_equalW"
else
    output=${output}"_varW"
fi
if [ "${shuffle_mode}" = "shuffle" ]; then
    output=${output}"_shuffle"
fi

# to plot distances
distance_output_file=${plot_dir}"/"${output}
if [ "${weights}" = "equal" ]; then
    "${plot_cmd[@]}" "${distance_output_file}" ${statDist} "${dssize2_list}" "${ceps2_list}" "${clap_n_distance_files} CLAP_EMPN" "${cgeom_n_distance_files} CGEOM_EMPN" "${cgauss_n_distance_files} CGAUSS_EMPN" "${ldp_gibu_distance_files} GIBU" 
else 
    "${plot_cmd[@]}" "${distance_output_file}" ${statDist} "${dssize2_list}" "${ceps2_list}" "${clap_n_distance_files} CLAP_EMPWN" "${cgeom_n_distance_files} CGEOM_EMPWN" "${cgauss_n_distance_files} CGAUSS_EMPWN" "${ldp_gibu_distance_files} GIBU" 
fi



echo " finished task" 

# use the following command to make the script executable ('a' for all users,'x' for executable ) 
# chmod a+x ./script.sh
