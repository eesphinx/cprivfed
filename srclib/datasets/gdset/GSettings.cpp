/*
 * GSettings.cpp
 *
 *  Created on: Jun 28, 2019
 *      Author: ehab
 */
#include <iostream>
#include <sstream>    // for reading from string as a file
#include <iomanip>
#include <ios>
#include "GSettings.h"

/** should be re-implemented by reading from a configuration file **/

std::string GSettings::get_setting_value(setting_name p_name) {
if (p_name == data_file)
	return "./gdata/Gowalla_totalCheckins.txt"; // data source file
    //return "./gdata/Brightkite_totalCheckins.txt"; // data source file
if (p_name == train_data_file)
	return "./gdata/gowalla-train.txt"; // train_data source file
if (p_name == test_data_file)
	return "./gdata/gowalla-test.txt"; // test_data source file
if (p_name == default_read_size)
	return "7000000";  // default number of records to be read from a dataset


// The following template is used when a flexible grid is created according to a
// specified cell size (cell_side) in the Settings.
if (p_name == map_script_flex_grid_template)
	return "./jscripts/map_script_flex_grid_template.html";

// used when a grid of fixed grid size (no. of rows & no. of columns) is created.
// In this case the cell size depends on the grid size.
if (p_name == map_script_fix_grid_template)
	return "./jscripts/fixed_grid/map_script_fixed_grid_template.html";

if (p_name == precision)
	return "10";

std::cerr << "error: unrecognized setting name";
exit(0);
return "";
}

/**
 * converts the given double value to a string with the given precision,
 * that is the number of decimal digits on the right of the fixed point.
 */
std::string GSettings::to_string_with_precision(double p_value, int p_precision) {
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(p_precision) << p_value;
	return oss.str();
}

/**
 * converts the given double value to a string with the default precision
 * in the settings.
 */
std::string GSettings::to_string_with_default_precision(double p_value) {
	return to_string_with_precision(p_value, std::stoi(GSettings::get_setting_value(precision)));
}

std::string GSettings::get_time() {
	auto current_time_point = std::chrono::system_clock::now();
	std::time_t current_time = std::chrono::system_clock::to_time_t(current_time_point);
	return std::ctime(&current_time);
	//if we need the elapsed seconds from start to end (both computed by std::chrono::system_clock::now())
	// std::chrono::duration<double> elapsed_seconds = end-start;
	// std::cout << elapsed_seconds.count();
}

