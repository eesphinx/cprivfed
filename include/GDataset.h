/*
 * GDataset.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_GDSET_GDATASET_H_
#define PRIV_PLAN_GDSET_GDATASET_H_


/*
A GDataset object is a list of check-in records that are read (according to boundaries and user criterion)
from a source file in Gowalla format. An instance is constructed with the following parameters:
data_file name:  the file containing the records in Gowalla format.
region boundaries: south, north, west, and east.
a user criteria: to obtain only records in which the user_id matches a given regular expression
read_size : number of records read from the data source file.

attributes:
1. region boundaries
2. default cell side: the region must form a grid of square cells having side def_cell_side
3. hsize,vsize: of the region grid
4. list of check-ins
5. number of users

methods:
1. write_to_file(p_data_file_name): write the data to a file in Gowalla format.
2. write_map_script: write a jscript file to plot the checkins in the dataset.
3. get_prior(cell_side): use cell_side to construct a grid and return a global prior on this grid.
	 Here cell_side must be multiples of def_cell_side, to get fine-grained prior.
4. write_all_priors_to_file(p_priors_file_name):

important: mapping between x,y coordinates of the cells and the lat/lng geographical coordinates is
necessary to construct priors.
x coordinate increases as we go in the east direction
y coordinate increases as we go towards the north.
*/

#include <vector>
#include <string>
#include <regex>
#include "GCheckin.h" // this is #included because we have a member vector<GCheckin>

class Prior;
class GridEmpData;
class UserStatistic;

class GDataset {
public:

    // if read_size is not specified it is set to a number larger than Gowalla dataset.
    GDataset(std::string p_data_file, double p_south, double p_north, double p_west, double p_east,
			double p_def_cell_side, std::regex p_user_id_regex, int p_read_size = 7000000);

    /**
     * get a vector containing the number of check-ins for every user having this number
     * above a given threshold. In this version, repeated points are considered one point.
     * */
    std::vector<UserStatistic> get_users_statistics_no_repeated(size_t p_user_checkins_threshold) const;


    /**
     * get a vector containing the number of check-ins for every user having this number
     * above the given threshold. In this version, repeated checkins for the same user
     * are considered different..
     * */
    std::vector<UserStatistic> get_users_statistics(int p_user_checkins_threshold) const;


	/** return an element in the dataset */
	GCheckin get_checkin(int index) const;

	/** returns the number of check-ins in the dataset. */
	int get_size() const;

	/** return the max distance in the region. */
	double get_max_distance() const;

	/** returns the number of users in the data set */
	int get_num_users() const;

	/** get information string of the data set */
	std::string get_info() const;

	/** write the data to the given file in Gwalla format. */
	void write_data_to_file(std::string p_file) const;


	/**
	 *  write all check-ins in the data set to a script that plots them on google maps
	 *  the file name should have the extension ".html". p_plot_type takes one of values {"MARKERS","HEATMAP"}
	 *  to determines whether to plot the markers of points, or a heatmap.
	 *  p_show_grid determines whether
	 *  or not to show the grip with the default cell_side on the map
	 **/
	void write_map_script(std::string p_script_file_name, std::string p_plot_type, bool p_show_grid) const;


	/**
	 * Return a prior constructed from the checkins. The dimensions of the prior are determined by p_cell_side.
	 **/
	Prior get_prior(double p_cell_side) const;

	/** construct a prior using the check-ins in the data set and the default cell_side */
	Prior get_prior() const;


	/** returns the global prior in a different way: constructs the prior of every user
	 * and then computes the global one by averaging users priors. */
	Prior get_average_prior(double p_cell_side) const;

	/** return the average prior using the check-ins in the data set in the above way and the default cell_side. */
	Prior get_average_prior() const;

	/**
	 * Return an EmpData object describing the original empirical data based on the checkins
	 * (without obfuscation). Each datum consists of a location and the frequency of of the occurrence
	 * of this location in the data. The granularity of the resulting empirical data is determined by
	 * p_cell_side (which must be multiples of def_cell_side).
	 **/
	GridEmpData get_grid_emp_data(double p_cell_side) const;


	/** construct a prior for every user in the data set and write it as a vector to the given file,
	 * i.e. in S format. parameters:
	 *  p_cell_side define the grid for which priors are constructed .
	 *  threshold: restricts the chosen users having at least this number of points.
	 **/
	void write_priors_of_users(std::string p_file, double p_cell_side, int p_user_checkins_threshold) const;

	/** write two priors for every user in p_ustat_list */
	void write_train_test_priors_of_users(const std::vector<UserStatistic>& p_ustat_list,
			std::string p_train_file, std::string p_test_file, double p_cell_side, double p_train_percentage) const;

	int get_number_of_distinct_checkins() const;

	/** return the hsize of the grid on the region based on the given p_cell_side */
	int get_hsize(double p_cell_side) const;

	/** return the vsize of the grid on the region based on the given p_cell_side */
	int get_vsize(double p_cell_side) const;

private:

	double south;
	double north;
	double west;
	double east;

	std::vector<GCheckin> checkins;

	double def_cell_side;
	// At the construction, the boundaries of the region is adjusted according to a
	// (default) def_cell_side to form a grid of square cells. def_cell_side is
	// also used to construct the grid on a the map when required to be shown.

    int num_users;  // number of users in the data set
    std::string dataset_info;

	/** read the data file sequentially and store the check-ins that are in the region and match p_user_id_regex */
	void read_checkins(std::string p_data_file, std::regex p_user_id_regex, int p_read_size);

	/** The following function adjusts the boundaries to be in form of grid where
	 * the side of every cell is def_cell_side */
	void adjust_boundaries();


	/** returns the horizontal index (int) of the cell (with length p_cell_side)
	 * containing the given point assuming that cell(0,0) is the most south-west cell */
	int get_cell_x_of_lng(double p_lng, double p_cell_side) const;

	/** returns the vertical index (int) of the cell (with length p_cell_side)
	 * containing the given point assuming that cell(0,0) is the most south-west cell */
	int get_cell_y_of_lat(double p_lat, double p_cell_side) const;

};




#endif /* PRIV_PLAN_GDSET_GDATASET_H_ */
