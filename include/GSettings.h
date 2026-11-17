/*
 * GSettings.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_GDSET_GSETTINGS_H_
#define PRIV_PLAN_GDSET_GSETTINGS_H_

#include <string>


// settings
enum setting_name { data_file, train_data_file, test_data_file, default_read_size,
	map_script_flex_grid_template, map_script_fix_grid_template, precision };


class GSettings {
public:

/** returns the value of the given setting in string format */
static std::string get_setting_value(setting_name p_name);

/**
 * converts the given double value to a string with the given precision,
 * that is the number of decimal digits on the right of the fixed point.
 */
static std::string to_string_with_precision(double p_value, int p_precision);

/**
 * converts the given double value to a string with the default precision
 * in the settings.
 */
static std::string to_string_with_default_precision(double p_value);

static std::string get_time();

};



#endif /* PRIV_PLAN_GDSET_GSETTINGS_H_ */
