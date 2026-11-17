# plots the estimation quality (statistical distance) for several protocols, 
# each one is specified by 8 files and a protocol label 
# Takes at least five arguments: the output file name, statistical distance, a list of leps2 values, 
# a list of ceps2 values, and the last one in the format "protocol_a_files label_a". 
# Additional 3 arguments for other protocols may be passed in the same format. 
#
# <output_file>      file name
# <statDist>         statistical distance
# "<leps2_list>"     list of leps2 values
# "<ceps2_list>"     list of ceps2 values resulting after shuffling
# "<protocol_a_files> labela"  protocol_a files and label
# "<protocol_b_files> labelb"  (optional) protocol_b files and label
# "<protocol_c_files> labelc"  (optional) protocol_c files and label
# "<protocol_d_files> labeld"  (optional) protocol_d files and label

if (ARGC < 5) {
    print "number of arguments = ".ARGC   
    print " Usage: gnuplot -c varying_leps2.gp <output_file> <statDist> ". \
          " ""<leps2_list>"" ""<ceps2_list>"" ". \
          " ""<protocol_a_files> labela"" ""<protocol_b_files> labelb (optional)"" ". \
          " ""<protocol_c_files> labelc (optional)"" ""<protocol_d_files> labeld (optional)"" "
    exit
}

output_file = ARG1
stat_distance = ARG2
leps2_list = ARG3
ceps2_list = ARG4
protocol_a = ARG5
protocol_b = (ARGC >= 6) ? ARG6 : ""
protocol_c = (ARGC >= 7) ? ARG7 : ""
protocol_d = (ARGC >= 8) ? ARG8 : ""

# Validation: 
# leps2_list must contain 8 elements: "<val1> <val2> ... <val8>"
if (words(leps2_list) != 8) {
    print "I got ", words(leps2_list) " values for the leps2. They must be exactly 8"
    exit
}

set output output_file.".pdf"

set terminal pdf linewidth 3 # use 3 for small figures, 2 for large figures
set tmargin 4.5   
set bmargin 4.5   # used if part of eps_g underneath the plot is hidden
set lmargin 10   # gives some space at the left for the y-label

set key font ",15" opaque # use additional "box" to put the key in a box
# set key right bottom
#set key center 
#set key spacing 1.5
# unset key

# y-settings
set logscale y
set format y "10^{%L}"
# set yrange [ 0.00001 : 0.1]  # adjust or comment out to display the full range or to unify the ranges of many plots
# y axis settings
if (stat_distance eq "SQD") {
    set ylabel font ",22" "squared error" #offset -1
    set ytics font ",15" # 1.0
}
if (stat_distance eq "TV") {
    set ylabel font ",22" "Total Variation" #offset -1
    set ytics font ",22" # 1.0
}
if (stat_distance eq "EMD") {
    set ylabel font ",22" "EMD" offset -1 # (km)" #offset 0
    set ytics font ",15"
}


# x-settings
set xlabel font ",22" "{/Symbol e}_2" offset 0,-0.5
set xtics border in scale 0,0 nomirror norotate # remove xticks
set xtics norangelimit
set xtics font ",22"
do for [i=1:words(leps2_list)] {
    set xtics add (word(leps2_list,i) i) # set xtics values
}


set linetype 17 lc rgb "#880000" #lw 2 pt 7   # variant A: dark red
set linetype 18 lc rgb "#FF7700" #lw 2 pt 5   # variant B: light red

# a function that returns the lineType that corresponds to a given protocol
lineType(p_protocol) = (p_protocol eq "CGAUSS_EMP") ? 1: \
                        (p_protocol eq "CGAUSS_EMPN") ? 1: \
                        (p_protocol eq "CGAUSS_EMPWN") ? 1: \
                        (p_protocol eq "CLAP_EMP") ? 2: \
                        (p_protocol eq "CLAP_EMPN") ? 2: \
                        (p_protocol eq "CLAP_EMPWN") ? 2: \
                        (p_protocol eq "CGEOM_EMP") ? 6: \
                        (p_protocol eq "CGEOM_EMPN") ? 6: \
                        (p_protocol eq "CGEOM_EMPWN") ? 6: \
                        (p_protocol eq "AVG_MECH_INV") ? 7: \
                        (p_protocol eq "AVG_MECH_INVP") ? 17: \
                        (p_protocol eq "AVG_MECH_INVN") ? 18: \
                        (p_protocol eq "GIBU") ? 8 :0
# 6: blue 7:red 5:yellow (unused) 

#Old: e.g. LAP+EMP+N. New: LAP+S&N
caption(p_protocol) = (p_protocol eq "CGAUSS_EMP") ?     "    GAUSS": \
                        (p_protocol eq "CGAUSS_EMPN") ?  "  GAUSS+N": \
                        (p_protocol eq "CGAUSS_EMPWN") ? " GAUSS+W+N": \
                        (p_protocol eq "CLAP_EMP") ?     "      LAP": \
                        (p_protocol eq "CLAP_EMPN") ?    "    LAP+N": \
                        (p_protocol eq "CLAP_EMPWN") ?   "   LAP+W+N": \
                        (p_protocol eq "CGEOM_EMP") ?    "     GEOM": \
                        (p_protocol eq "CGEOM_EMPN") ?   "   GEOM+N": \
                        (p_protocol eq "CGEOM_EMPWN") ?  "  GEOM+W+N": \
                        (p_protocol eq "AVG_MECH_INV") ? "ObfAndCount+INV": \
                        (p_protocol eq "AVG_MECH_INVP") ? "ObfAndCount+INV+P": \
                        (p_protocol eq "AVG_MECH_INVN") ? "ObfAndCount+INV+N": \
                        (p_protocol eq "GIBU") ?         "   kRR+GIBU" : "Unknown"

# show the top labels if ceps2_list is nonempty
unset label
if (words(ceps2_list)>0) {
# set "ceps" (label 9) at the top and list the values underneath
set label 9 sprintf("{/Symbol e}^c_2") \
    at 4.5, graph 1.15 \
    center \
    font ",18" \
    tc rgb "black" \
    front \
    offset 0, 0
do for [i=1:8] {
    x_pos = i  # Group position on x-axis
    # cmd = sprintf("./krr_clones_feldman2020 %g %g %g %g", real(leps1), real(delta), real(shuffle_dssize), 100)
    # eps1 = (shuff eq "shuffle") ? real(system(cmd)) : real(leps1)
    # cmd = sprintf("./krr_clones_feldman2020 %g %g %g %g", real(word(leps2_range,i)), real(delta), real(shuffle_dssize), real(alphabet_size))
    # eps2 = (shuff eq "shuffle") ? real(system(cmd)) : real(word(leps2_range,i))
    # Position label above the boxplot group
    set label i sprintf("%.3f", real(word(ceps2_list,i))) \
        at x_pos, graph 1.05 \
        center \
        font ",15" \
        tc rgb "black" \
        front \
        offset 0, 0
}
}


# construct the plot command dynamically
plot_cmd = "plot "

if (protocol_a ne "") {
    if (words(protocol_a) != 9) {
        print "invalid first protocol data. Must have 8 file names"
        exit
    }
    protocol_a_label = word(protocol_a,9)
    data_protocol_a_01 = word(protocol_a,1)
    data_protocol_a_02 = word(protocol_a,2)
    data_protocol_a_03 = word(protocol_a,3)
    data_protocol_a_04 = word(protocol_a,4)
    data_protocol_a_05 = word(protocol_a,5)
    data_protocol_a_06 = word(protocol_a,6)
    data_protocol_a_07 = word(protocol_a,7)
    data_protocol_a_08 = word(protocol_a,8)

    plot_cmd = plot_cmd. \
    "data_protocol_a_01 using (1-s):1 lt ".lineType(protocol_a_label)." title "."\"".caption(protocol_a_label)."\"".", ". \
    "data_protocol_a_02 using (2-s):1 lt ".lineType(protocol_a_label)." notitle, ". \
    "data_protocol_a_03 using (3-s):1 lt ".lineType(protocol_a_label)." notitle, ". \
    "data_protocol_a_04 using (4-s):1 lt ".lineType(protocol_a_label)." notitle, ". \
    "data_protocol_a_05 using (5-s):1 lt ".lineType(protocol_a_label)." notitle, ". \
    "data_protocol_a_06 using (6-s):1 lt ".lineType(protocol_a_label)." notitle, ". \
    "data_protocol_a_07 using (7-s):1 lt ".lineType(protocol_a_label)." notitle, ". \
    "data_protocol_a_08 using (8-s):1 lt ".lineType(protocol_a_label)." notitle,"
}

if (protocol_b ne "") {
    if (words(protocol_b) != 9) {
        print "invalid protocol_b data. Must have 8 file names\n"
        exit
    }
    protocol_b_label = word(protocol_b,9)
    data_protocol_b_01 = word(protocol_b,1)
    data_protocol_b_02 = word(protocol_b,2)
    data_protocol_b_03 = word(protocol_b,3)
    data_protocol_b_04 = word(protocol_b,4)
    data_protocol_b_05 = word(protocol_b,5)
    data_protocol_b_06 = word(protocol_b,6)
    data_protocol_b_07 = word(protocol_b,7)
    data_protocol_b_08 = word(protocol_b,8)

    plot_cmd = plot_cmd. \
    "data_protocol_b_01 using (1):1 lt ".lineType(protocol_b_label)." title "."\"".caption(protocol_b_label)."\"".", ". \
    "data_protocol_b_02 using (2):1 lt ".lineType(protocol_b_label)." notitle, ". \
    "data_protocol_b_03 using (3):1 lt ".lineType(protocol_b_label)." notitle, ". \
    "data_protocol_b_04 using (4):1 lt ".lineType(protocol_b_label)." notitle, ". \
    "data_protocol_b_05 using (5):1 lt ".lineType(protocol_b_label)." notitle, ". \
    "data_protocol_b_06 using (6):1 lt ".lineType(protocol_b_label)." notitle, ". \
    "data_protocol_b_07 using (7):1 lt ".lineType(protocol_b_label)." notitle, ". \
    "data_protocol_b_08 using (8):1 lt ".lineType(protocol_b_label)." notitle,"
}

if (protocol_c ne "") {
    if (words(protocol_c) != 9) {
        print "invalid protocol_c data. Must have 8 file names"
        exit
    }
    protocol_c_label = word(protocol_c,9)
    data_protocol_c_01 = word(protocol_c,1)
    data_protocol_c_02 = word(protocol_c,2)
    data_protocol_c_03 = word(protocol_c,3)
    data_protocol_c_04 = word(protocol_c,4)
    data_protocol_c_05 = word(protocol_c,5)
    data_protocol_c_06 = word(protocol_c,6)
    data_protocol_c_07 = word(protocol_c,7)
    data_protocol_c_08 = word(protocol_c,8)

    plot_cmd = plot_cmd. \
    "data_protocol_c_01 using (1+s):1 lt ".lineType(protocol_c_label)." title "."\"".caption(protocol_c_label)."\"".", ". \
    "data_protocol_c_02 using (2+s):1 lt ".lineType(protocol_c_label)." notitle, ". \
    "data_protocol_c_03 using (3+s):1 lt ".lineType(protocol_c_label)." notitle, ". \
    "data_protocol_c_04 using (4+s):1 lt ".lineType(protocol_c_label)." notitle, ". \
    "data_protocol_c_05 using (5+s):1 lt ".lineType(protocol_c_label)." notitle, ". \
    "data_protocol_c_06 using (6+s):1 lt ".lineType(protocol_c_label)." notitle, ". \
    "data_protocol_c_07 using (7+s):1 lt ".lineType(protocol_c_label)." notitle, ". \
    "data_protocol_c_08 using (8+s):1 lt ".lineType(protocol_c_label)." notitle,"
}

if (protocol_d ne "") {
    if (words(protocol_d) != 9) {
        print "invalid protocol_d data. Must have 8 file names"
        exit
    }
    protocol_d_label = word(protocol_d,9)
    data_protocol_d_01 = word(protocol_d,1)
    data_protocol_d_02 = word(protocol_d,2)
    data_protocol_d_03 = word(protocol_d,3)
    data_protocol_d_04 = word(protocol_d,4)
    data_protocol_d_05 = word(protocol_d,5)
    data_protocol_d_06 = word(protocol_d,6)
    data_protocol_d_07 = word(protocol_d,7)
    data_protocol_d_08 = word(protocol_d,8)

    plot_cmd = plot_cmd. \
    "data_protocol_d_01 using (1):1 lt ".lineType(protocol_d_label)." title "."\"".caption(protocol_d_label)."\"".", ". \
    "data_protocol_d_02 using (2):1 lt ".lineType(protocol_d_label)." notitle, ". \
    "data_protocol_d_03 using (3):1 lt ".lineType(protocol_d_label)." notitle, ". \
    "data_protocol_d_04 using (4):1 lt ".lineType(protocol_d_label)." notitle, ". \
    "data_protocol_d_05 using (5):1 lt ".lineType(protocol_d_label)." notitle, ". \
    "data_protocol_d_06 using (6):1 lt ".lineType(protocol_d_label)." notitle, ". \
    "data_protocol_d_07 using (7):1 lt ".lineType(protocol_d_label)." notitle, ". \
    "data_protocol_d_08 using (8):1 lt ".lineType(protocol_d_label)." notitle,"
}


# boxplot style
set style boxplot nooutliers
set style data boxplot
set boxwidth 0.3 absolute
set style fill solid 0.25 border # lt -1


# this adds vertical separator lines at points 1.5, 2.5, etc
set arrow from 1.5, graph 0 to 1.5, graph 1 nohead dt "."
set arrow from 2.5, graph 0 to 2.5, graph 1 nohead dt "."
set arrow from 3.5, graph 0 to 3.5, graph 1 nohead dt "."
set arrow from 4.5, graph 0 to 4.5, graph 1 nohead dt "."
set arrow from 5.5, graph 0 to 5.5, graph 1 nohead dt "."
set arrow from 6.5, graph 0 to 6.5, graph 1 nohead dt "."
set arrow from 7.5, graph 0 to 7.5, graph 1 nohead dt "."

set grid ytics lt 0
set xrange [ 0.5 : 8.5]  # to display the full range

s=0.2 # gap between neighbor boxes 

str_cmd = substr(plot_cmd, 1, strlen(plot_cmd)-1) # remove the last ","
eval(plot_cmd)


# reuse from the following if needed
#set key box opaque font "Verdana,15"
#set title font "Verdana,15" "using xxx mechanism for obfuscating \n user data in a planar region (24 cells west-east, and 12 #cells south-north)"
#set xlabel font "Verdana,15" "epsilon"
#set xlabel font ",20" "{/Symbol e}_g"
#set ylabel font "Verdana,15" "Earth Mover"s distance"
#set xtics font "Verdana,15"
#set ytics font "Verdana,15"
#set title "using KRR mechanism for obfuscating \n user data in a planar region (24 cells west-east, and 12 cells south-north)"

# y axis settings
# uncomment for TV:
#set ylabel font ",20" "Total Variation" offset 0
#set yrange [ 0 : 1]  # to display the full range
#set ytics 0.1
# uncomment for EMD:
#set ylabel font ",20" ytitle offset 0
#set yrange [ 0 : 25]  # to display the full range
#set ytics 1



