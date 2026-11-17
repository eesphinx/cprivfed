
#include <regex>
#include <vector>
#include <algorithm>
#include <sstream>    // for reading from string as a file
#include <fstream>    // for reading and writing to files
#include <iostream>   // for cerr, cout
#include <cmath>
#include "GDataset.h"
#include "UserStatistic.h"
#include "GCheckin.h"
#include "GSettings.h"
#include "Distance.h"
#include "Prior.h"
#include "GridEmpData.h"
#include "GridObservation.h"
#include "Mechanism.h"


/** if read_size is not specified it is set to a number larger than Gowalla dataset */
GDataset::GDataset(std::string p_data_file, double p_south, double p_north, double p_west, double p_east,
		double p_def_cell_side, std::regex p_user_id_regex, int p_read_size)
:south(p_south), north(p_north), west(p_west), east(p_east), def_cell_side(p_def_cell_side) {
	adjust_boundaries();
	num_users = 0;
	read_checkins(p_data_file,p_user_id_regex,p_read_size);
}

/**
 * read the data file sequentially and store the check-ins that are in the
 * region and match p_user_id_regex
 */
void GDataset::read_checkins(std::string p_data_file, std::regex p_user_id_regex, int p_read_size) {
	int read_row_count=0;     // number of valid rows read from in_data_file
	std::ifstream pfile(p_data_file);
	if (pfile.is_open())
	{
		// process the contents of the file
		std::string line;
		// summary data
		std::string old_user_id = "*"; // holds the previous user in the output to be compared with the new (curr) user.  This is important to count the users in the dataset
		while (getline(pfile,line)) {
			if (read_row_count >= p_read_size)
				break;
			if (line.length()==0 || line.front()=='#')
				continue;
			read_row_count++;
			GCheckin chkin(line); // construct a checkin object from a line
			//cout << line << endl;
			//cout << chkin.user_id << ' ' << chkin.lat << ' ' << chkin.lng << endl;
			if (chkin.lat.length()==0 || chkin.lng.length()==0)
				continue;
			if (stod(chkin.lat) < south || stod(chkin.lat) > north
					|| stod(chkin.lng) < west || stod(chkin.lng) > east)
				continue;
			if (!regex_match(chkin.user_id, p_user_id_regex)) // handle regular expressions
				continue;
			// here we have a valid checkin
			checkins.push_back(chkin);
			if (chkin.user_id.compare(old_user_id)!=0) {
				old_user_id = chkin.user_id;
				num_users++;
			}
		}
		pfile.close();
	}
	else {
		std::cerr << "Unable to open dataset file " << p_data_file;
		exit(0);
	}

	std::ostringstream strs;

	strs << "# Boundaries: \n";
	strs << "# South: "<< south <<             //  String.format("%03.4f", south)+
			",   North: "<<  north << "\n";                 //String.format("%03.4f", north)+"\n";
	strs << "# West:  "<< west <<              // String.format("%03.4f", west)+
			",    East:  "<< east << "\n";                  // String.format("%03.4f", east)+"\n";
	strs << "# Number of checkins = " << checkins.size() <<"\n";
	strs << "# filtered from  "<< read_row_count<< "  input rows"<<"\n";
	strs << "# Total number of users in the region = " << num_users;
	dataset_info = strs.str();
}



/**
 * get a vector containing the number of check-ins for every user having this number
 * above the given threshold. In this version, repeated points are considered one point.
 * */
std::vector<UserStatistic> GDataset::get_users_statistics_no_repeated(size_t p_user_checkins_threshold) const{
	std::vector<UserStatistic> retList;  // return vector
	UserStatistic us;
	std::string old_user_id = checkins[0].user_id;   // (*checkins_iterator).user_id;

	std::vector<std::string> old_user_points;  // points of the same user
	old_user_points.push_back(checkins[0].lat + " "+ checkins[0].lat);

	for (size_t k=1;k<checkins.size();k++) {
		// check if a new user is the same as before
		if ( checkins[k].user_id== old_user_id) {
			// if the point is already encountered for this user
			if (find( old_user_points.begin(), old_user_points.end(),
					checkins[k].lat + " "+ checkins[k].lat ) == old_user_points.end()) {
				old_user_points.push_back(checkins[k].lat + " "+ checkins[k].lat);
			}
		}
		// here we have a new user
		else {
			// process the info of the old user first (if bigger than threshold)
			if (old_user_points.size() >= p_user_checkins_threshold) {
				us.set(old_user_id, old_user_points.size());
				retList.push_back(us);
			}
			// process the new user
			old_user_id = checkins[k].user_id;
			old_user_points.clear();
			old_user_points.push_back(checkins[k].lat + " "+ checkins[k].lat);
		}
	}
	// process the info of the last user (if bigger than threshold)
	if (old_user_points.size() >= p_user_checkins_threshold) {
		us.set(old_user_id,old_user_points.size());
		retList.push_back(us);
	}
	return retList;
}

/**
 * get a vector containing the number of check-ins for every user having this number
 * above the given threshold. In this version, repeated checkins for the same user
 * are considered different..
 * Although this function does not change the state of the object we do not specify it as [const]
 * because we access an iterator on the checkins.
 * */
std::vector<UserStatistic> GDataset::get_users_statistics(int p_user_checkins_threshold) const {

	std::vector<UserStatistic> retList;  // return vector
	UserStatistic us;
	std::string old_user_id = checkins[0].user_id;
	int old_user_checkins=1;   // note here we don't need to keep a list of checkins of the old user, just track the count

	for (size_t k=1;k<checkins.size();k++) {
		// check if a new user is the same as before
		if ( checkins[k].user_id == old_user_id)
			old_user_checkins++;
		else {
			// process the info of the old user first (if bigger than threshold)
			if (old_user_checkins >= p_user_checkins_threshold) {
				us.set(old_user_id, old_user_checkins);
				retList.push_back(us);
			}
			// process the new user
			old_user_id = checkins[k].user_id;
			old_user_checkins=1;
		}
	}
	// process the info of the old user
	if (old_user_checkins >= p_user_checkins_threshold) {
		us.set(old_user_id, old_user_checkins);
		retList.push_back(us);
	}
	return retList;
}



/** return an element in the dataset */
GCheckin GDataset::get_checkin(int index) const {
	return checkins.at(index);
}

/** returns the number of check-ins in the dataset. */
int GDataset::get_size() const {
	return checkins.size();
}

/** return the max distance in the region. */
double GDataset::get_max_distance() const {
	return Distance::get_distance_on_earth(north, east, south, west);
}

/** returns the number of users in the data set */
int GDataset::get_num_users() const {
	return num_users;
}

/** get information string of the data set */
std::string GDataset::get_info() const {
	return dataset_info;
}


/** write the data to the given file in Gwalla format */
void GDataset::write_data_to_file(std::string p_file) const {
	std::ofstream pfile (p_file);
	if (pfile.is_open()) {
		// write data
		for (size_t k=0; k< checkins.size(); k++)
			pfile << checkins[k].user_id << "\t" << checkins[k].datetime+"\t" << checkins[k].lat << "\t" << checkins[k].lng << "\t" << checkins[k].loc_id << "\n";
		pfile << "\n\n"+ dataset_info;
		pfile.close();

	}
	else {
		std::cerr << "Unable to open file for writing the dataset";
		exit(0);
	}
}

// write all check-ins in the data set to a script that plots them on google maps
// the file name should have the extension ".html". p_plot_type takes one of values {"MARKERS","HEATMAP"}
// to determines whether to plot the markers of points, or a heatmap.
//p_show_grid determines whether
// or not to show the grip with the default cell_side on the map.
void GDataset::write_map_script(std::string p_script_file_name, std::string p_plot_type, bool p_show_grid) const {
	std::string map_script_flex_grid_template = "./jscripts/map_script_flex_grid_template.html";
	std::string map_script_fix_grid_template = "./jscripts/fixed_grid/map_script_fixed_grid_template.html";

	std::ifstream in_file(map_script_flex_grid_template);
	if (!in_file.is_open()) {
		std::cerr << "Unable to open script template file: " << p_script_file_name;
		exit(0);
	}
	std::ofstream out_file(p_script_file_name);
	if (!out_file.is_open()) {
		std::cerr << "Unable to open file for writing the dataset ";
		exit(0);
	}

	std::string line;
	// Here we copy the template_script to out_script until "// data" is found
	std::string replacement_tag = "// to be replaced";
	while (getline(in_file,line) &&
			line.substr(0, replacement_tag.length()) != replacement_tag)
		out_file << line << "\n";

	// Here we start writing the out_data
	out_file << "var rect_boundaries = {\n"
			<< "      south: " << south << ",\n"
			<< "      north: " << north << ",\n"
			<< "      west : " << west  << ",\n"
			<< "      east : " << east  << ",\n"
			<< "    };\n"
			<< "//\n"
			<< "// checkins\n"
			<< "var locs = [\n";


	for (size_t k=0; k< checkins.size(); k++)
		out_file << "            {lat:" << checkins[k].lat << " , lng:" << checkins[k].lng << "},\n";

	//  At this point data has been filtered, now write the variables.
	out_file << "];\n"
			<< "// var loc_count =" << checkins.size() << "; \n"
			<< "//\n"
			<< "// grid settings: \n"
			<< "var grid_required = " << p_show_grid <<"; \n"
			<< "var plot_type = \"" << p_plot_type  << "\"; \n"
			<< "var d_u =" << def_cell_side << "; \n";

	// Continue to copy from the template script to the out_script
	while (getline(in_file,line)) out_file << line << "\n";
	out_file << std::endl;   // flush the buffer

	in_file.close();
	out_file.close();
}

/**
 * Return a prior constructed from the checkins. The dimensions of the prior are determined by p_cell_side.
 **/
Prior GDataset::get_prior(double p_cell_side) const {
	double scaling_ratio = def_cell_side/p_cell_side;
	if (scaling_ratio != (int) scaling_ratio) {
		std::cerr << "unable to get a prior from dataset:\n";
		std::cerr << "The default cell side must be multiples of the given cell side";
		exit(0);
	}

	int hsize = get_hsize(p_cell_side);
	int vsize = get_vsize(p_cell_side);

	// create a row vector for prior (elements are stored in linear order) and initialize to 0
	DoubleRow prior_v(hsize*vsize); prior_v.fillWith(0);
	for (size_t k=0; k< checkins.size(); k++) {
		int cell_x = get_cell_x_of_lng(   stod(checkins[k].lng),  p_cell_side);
		int cell_y = get_cell_y_of_lat(   stod(checkins[k].lat),  p_cell_side);
		prior_v.put(cell_y*hsize + cell_x, prior_v.get(cell_y*hsize + cell_x) +1 );
	}
	// normalize
	prior_v.divi(checkins.size());
	Prior pr(prior_v, hsize, vsize);

	pr.add_info(dataset_info);
	return pr;
}

// return a prior using the check-ins in the data set and
// the default cell_side
Prior GDataset::get_prior() const {
	return get_prior(def_cell_side);
}


// returns the global prior in a different way: constructs the prior of every
// user and then computes the global one by averaging users priors.
Prior GDataset::get_average_prior(double p_cell_side) const {
	double scaling_ratio = def_cell_side/p_cell_side;
	if (scaling_ratio != (int) scaling_ratio) {
		std::cerr << "unable to construct average prior from dataset:\n";
		std::cerr << "The default cell side must be multiples of the given cell side";
		exit(0);
	}

	int hsize = get_hsize(p_cell_side); int vsize = get_vsize(p_cell_side);
	// create a row vectors for the global and current user priors
	DoubleRow gprior_v(hsize*vsize); gprior_v.fillWith(0);
	DoubleRow uprior_v(hsize*vsize); uprior_v.fillWith(0);

	// summary data
	int out_users_count=0;  // count the number of users in the output

	// handle initial user
	int cell_x = get_cell_x_of_lng(   stod(checkins[0].lng),  p_cell_side);
	int cell_y = get_cell_y_of_lat(   stod(checkins[0].lat),  p_cell_side);
	std::string old_user_id = checkins[0].user_id;
	int old_user_checkins=1;
	// update the user's prior matrix
	uprior_v.put(cell_y*hsize+cell_x, uprior_v.get(cell_y*hsize+cell_x)+1);

	for (size_t k = 1; k < checkins.size(); k++) {
		// determine the cell of the curr_coordinates.
		int cell_x = get_cell_x_of_lng(   stod(checkins[k].lng),  p_cell_side);
		int cell_y = get_cell_y_of_lat(   stod(checkins[k].lat),  p_cell_side);

		// check if a new user is the same as before
		if ( checkins[k].user_id == old_user_id) {
			old_user_checkins++;
			uprior_v.put(cell_y*hsize+cell_x, uprior_v.get(cell_y*hsize+cell_x)+1);
		}
		else {
			// process the info of the old user
			out_users_count++;
			// normalize the user's vector and update the global prior vector
			uprior_v.divi(old_user_checkins);
			gprior_v.addi(uprior_v);

			// now process new user
			// reset the uprior array
			uprior_v.fillWith(0);
			old_user_id = checkins[k].user_id;
			old_user_checkins=1;
			// update uprior
			uprior_v.put(cell_y*hsize+cell_x, uprior_v.get(cell_y*hsize+cell_x)+1);
		}
	}
	// process the last user
	out_users_count++;
	// normalize the user's vector and update the global prior vector
	uprior_v.divi(old_user_checkins);
	gprior_v.addi(uprior_v);

	// average the global prior vector entries (over users), and create a prior using gprior_v
	gprior_v.divi(out_users_count);
	Prior gprior(gprior_v, hsize, vsize);
	gprior.add_info(dataset_info);

	return gprior;
}

/** return the average prior using the check-ins in the data set in the above way and the default cell_side. */
Prior GDataset::get_average_prior() const {
	return get_average_prior(def_cell_side);
}


/**
 * Return a GridEmpData object describing the original grid empirical data based on the checkins
 * (without obfuscation). Each datum consists of a location and the frequency of of the occurrence
 * of this location in the data. The granularity of the resulting empirical data is determined by
 * p_cell_side (which must be multiples of def_cell_side).
 **/
GridEmpData GDataset::get_grid_emp_data(double p_cell_side) const {
	double scaling_ratio = def_cell_side/p_cell_side;  // allows finer grids but not coarser.
	if (scaling_ratio != (int) scaling_ratio) {
		std::cerr << "unable to construct empirical data from dataset:\n";
		std::cerr << "The default cell side must be multiples of the given cell side";
		exit(0);
	}

	int hsize = get_hsize(p_cell_side);
	int vsize = get_vsize(p_cell_side);

	// create an array to record the frequency of every cell. Cells are ordered in the array
	// in the linear order.
	std::vector<int> freq_array(hsize*vsize , 0); // create array and initialize to 0.  
	for (size_t k=0; k< checkins.size(); k++) {
		//update the prior_matrix;
		int cell_x = get_cell_x_of_lng(   stod(checkins[k].lng),  p_cell_side);
		int cell_y = get_cell_y_of_lat(   stod(checkins[k].lat),  p_cell_side);
		freq_array[cell_y * hsize + cell_x] = freq_array[cell_y * hsize + cell_x]+1;
	}
	GridEmpData grid_data(freq_array, hsize, vsize);
	return grid_data;
}


/**
 * construct a prior for every user in the data set and write it as a vector to
 * the given file, i.e. in S format.
 * parameters:
 * p_cell_side define the grid for which priors are constructed .
 * threshold: restricts the chosen users having at least this number of points.
 */
void GDataset::write_priors_of_users(std::string p_file, double p_cell_side, int p_user_checkins_threshold) const {
	std::ofstream pfile(p_file);
	if (!pfile.is_open()) {
		std::cerr << "Unable to open file for writing users priors";
		exit(0);
	}

	int hsize = get_hsize(p_cell_side);
	int vsize = get_vsize(p_cell_side);

	DoubleRow uprior_v(hsize*vsize); uprior_v.fillWith(0);

	// summary data
	int out_users_count=0;  // count the number of users in the output

	// process the first check-in
	// handle initial user
	int cell_x = get_cell_x_of_lng(   stod(checkins[0].lng),  p_cell_side);
	int cell_y = get_cell_y_of_lat(   stod(checkins[0].lat),  p_cell_side);
	std::string old_user_id = checkins[0].user_id;
	int old_user_checkins=1;
	// update the user's prior array: add the new position
	uprior_v.put(cell_y*hsize+cell_x, uprior_v.get(cell_y*hsize+cell_x)+1);

	for (size_t k = 1 ; k<checkins.size(); k++ ) {
		// determine the cell of the curr_coordinates.
		int cell_x = get_cell_x_of_lng(   stod(checkins[k].lng),  p_cell_side);
		int cell_y = get_cell_y_of_lat(   stod(checkins[k].lat),  p_cell_side);
		// check if a new user is the same as before
		if (checkins[k].user_id == old_user_id) {
			old_user_checkins++;
			uprior_v.put(cell_y*hsize+cell_x, uprior_v.get(cell_y*hsize+cell_x)+1);
		}
		else {
			// process the info of the old user
			if (old_user_checkins >= p_user_checkins_threshold) {
				out_users_count++;
				// here we normalize the user's prior vector, create the prior, and write the prior to file
				uprior_v.divi(old_user_checkins);
				Prior uprior(uprior_v, hsize, vsize);
				pfile << "# user " << old_user_id << ": no of checkins = " << old_user_checkins << "\n";
				pfile << uprior.toString() << "\n";
			}
			// now we have a new user: reset the user prior vector
			uprior_v.fillWith(0);
			old_user_id = checkins[k].user_id;
			old_user_checkins=1;
			// add the new position
			uprior_v.put(cell_y*hsize+cell_x, uprior_v.get(cell_y*hsize+cell_x)+1);
		}
	}
	// process the last user
	if (old_user_checkins >= p_user_checkins_threshold) {
		out_users_count++;
		//  here we normalize the user's matrix
		uprior_v.divi(old_user_checkins);
		// create a prior and write its data to the file
		Prior uprior(uprior_v, hsize,vsize);
		pfile << "# user " << old_user_id  << ": no of checkins = " << old_user_checkins << "\n";
		pfile << uprior.toString() << "\n";
		// reset the user prior vector
		uprior_v.fillWith(0);
	}
	// prepare the tail summary data to be written to out_data_file
	pfile << "\n\n" << dataset_info << "\n" <<
			"# Number of users (with at least " << p_user_checkins_threshold
			<< " checkins) = " << out_users_count<< std::endl;
}



/** write two priors for every user in p_ustat_list */
void GDataset::write_train_test_priors_of_users(const std::vector<UserStatistic>& p_ustat_list,
		std::string p_train_file, std::string p_test_file, double p_cell_side, double p_train_percentage) const {

	std::ofstream train_file(p_train_file);
	if (!train_file.is_open()) {
		std::cerr << "Unable to open file for writing users training priors ";
		exit(0);
	}
	std::ofstream test_file(p_test_file);
	if (!test_file.is_open()) {
		std::cerr << "Unable to open file for writing users testing priors ";
		exit(0);
	}

	int hsize = get_hsize(p_cell_side);
	int vsize = get_vsize(p_cell_side);

	// create vectors for the the current user prior
	DoubleRow uprior_train_v(hsize*vsize); uprior_train_v.fillWith(0);
	DoubleRow uprior_test_v(hsize*vsize); uprior_test_v.fillWith(0);

	// summary data
	//int out_users_count=0;  // count the number of users in the output

	std::string old_user_id ="*";
	// holds the previous user in the output to be compared with the new (curr) user.
	int old_user_train_checkins =0;
	int old_user_test_checkins =0;


	// check if the user is in the list
	bool old_user_in_list = false;
	size_t old_user_index_in_list = 0;
	while (old_user_index_in_list < p_ustat_list.size() && p_ustat_list[old_user_index_in_list].user_id != old_user_id)
		old_user_index_in_list++;
	if (old_user_index_in_list<  p_ustat_list.size())
		old_user_in_list = true;



	// counting the number of check-ins for the old user.
	for (size_t it = 0; it < checkins.size(); it++) {
		// determine the cell of the curr_coordinates.
		int cell_x = get_cell_x_of_lng(   stod(checkins[it].lng),  p_cell_side);
		int cell_y = get_cell_y_of_lat(   stod(checkins[it].lat),  p_cell_side);
		// check if a new user is the same as before
		if (checkins[it].user_id==old_user_id) {
			//update the user's prior_vector;
			if (old_user_in_list) {
				if (old_user_train_checkins+old_user_test_checkins+1 < p_train_percentage* p_ustat_list[old_user_index_in_list].count) {
					uprior_train_v.put(cell_y*hsize+cell_x, uprior_train_v.get(cell_y*hsize+cell_x)+1);
					old_user_train_checkins++;
				}
				else {
					uprior_test_v.put(cell_y*hsize+cell_x, uprior_test_v.get(cell_y*hsize+cell_x)+1);
					old_user_test_checkins++;
				}
			}
		}
		else {
			// process the info of the old user
			if (old_user_in_list) {
				// here we normalize the user's train and test vectors, create priors, and write data to file
				uprior_train_v.divi(old_user_train_checkins);
				uprior_test_v.divi(old_user_test_checkins);
				Prior train_uprior(uprior_train_v, hsize,vsize);
				Prior test_uprior(uprior_test_v, hsize,vsize);
				train_file << "# user "+ old_user_id << ": no of checkins = " << old_user_train_checkins << "\n";
				train_file << train_uprior.toString() << "\n";
				test_file << "# user " << old_user_id  << ": no of checkins = " << old_user_test_checkins << "\n";
				test_file << test_uprior.toString() << "\n";
			}
			// now we have a new user: reset the vectors
			uprior_train_v.fillWith(0);
			uprior_test_v.fillWith(0);
			old_user_id = checkins[it].user_id;

			// initialize data for the new user
			old_user_train_checkins =0;
			old_user_test_checkins =0;
			// check if the user is in the list

			// check if the user is in the list
			bool old_user_in_list = false;
			size_t old_user_index_in_list = 0;
			while (old_user_index_in_list < p_ustat_list.size() && p_ustat_list[old_user_index_in_list].user_id != old_user_id)
				old_user_index_in_list++;
			if (old_user_index_in_list < p_ustat_list.size())
				old_user_in_list = true;

			old_user_train_checkins=1;
			if (old_user_in_list)
				uprior_train_v.put(cell_y*hsize+cell_x, uprior_train_v.get(cell_y*hsize+cell_x)+1);
		}
	}
	// to display the information of the last user
	// process the info of the old user
	if (old_user_in_list) {
		// normalize the user's prior vectors, create and write data to file
		uprior_train_v.fillWith(0);
		uprior_test_v.fillWith(0);
		Prior train_uprior(uprior_train_v, hsize,vsize);
		Prior test_uprior(uprior_test_v, hsize, vsize);
		train_file << "# user " << old_user_id  << ": no of checkins = " << old_user_train_checkins << "\n";
		train_file << train_uprior.toString() << "\n";
		test_file << "# user "<< old_user_id <<": no of checkins = " << old_user_test_checkins << "\n";
		test_file << test_uprior.toString() << "\n";
	}
	// now we have a new user : reset the vectors
	uprior_train_v.fillWith(0);
	uprior_test_v.fillWith(0);
}



/**
 * The following function adjusts the boundaries to be in form of grid
 * where the side of every cell is def_cell_side
 */
void GDataset::adjust_boundaries() {
	double d_ns = Distance::get_distance_on_earth(north, east, south, east);
	double d_ew = Distance::get_distance_on_earth(north, east, north, west);
	int vsize = (int) round(d_ns/def_cell_side);
	int hsize = (int) round(d_ew/def_cell_side);
	north = south+ (north-south)*(def_cell_side/d_ns)*vsize;
	east = west+ (east-west)*(def_cell_side/d_ew)*hsize;
}

int GDataset::get_hsize(double p_cell_side) const {
	double d_ew = Distance::get_distance_on_earth(north, east, north, west);
	return (int) round(d_ew/p_cell_side);
}
int GDataset::get_vsize(double p_cell_side) const {
	double d_ns = Distance::get_distance_on_earth(north, east, south, east);
	return (int) round(d_ns/p_cell_side);
}


/** return the number of distinct points in the dataset */
int GDataset::get_number_of_distinct_checkins() const {
	std::vector<std::string> points;  // points
	for (size_t it= 0; it < checkins.size(); it++)
		// if the point is already encountered for this user, i.e. not found in points
		if (find(points.begin(), points.end(),checkins[it].lat +" "+ checkins[it].lng ) == points.end())
			points.push_back(checkins[it].lat + " "+ checkins[it].lng);
	return points.size();
}


/**
 * returns the horizontal index (int) of the cell (with length p_cell_side)
 * containing the given point assuming that cell(0,0) is the most south-west cell.
 */
int GDataset::get_cell_x_of_lng(double p_lng, double p_cell_side) const {
	int hsize = get_hsize(p_cell_side);
	double h_spacing = (east - west)/hsize;
	return (int)((p_lng- west)/h_spacing);
}

/**
 * returns the vertical index (int) of the cell (with length p_cell_side)
 * containing the given point assuming that cell(0,0) is the most south-west cell
 */
int GDataset::get_cell_y_of_lat(double p_lat, double p_cell_side) const {
	int vsize = get_vsize(p_cell_side);
	double v_spacing = (north - south)/vsize;
	return (int)((p_lat- south)/v_spacing);
}


