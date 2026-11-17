/*
 * GCheckin.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_GDSET_GCHECKIN_H_
#define PRIV_PLAN_GDSET_GCHECKIN_H_

#include <string>

class GCheckin {
public:
	std::string user_id;
	std::string datetime;
	std::string lat;
	std::string lng;
	std::string loc_id;

	// construct a checkin record from a (string) line
	GCheckin(const std::string p_data_line);

	//GCheckin();

	//GCheckin(std::string p_user_id, std::string p_datetime, std::string p_lat, std::string p_lng, std::string p_loc_id);

	//void copy_data_from_line(const std::string p_data_line);
};




#endif /* PRIV_PLAN_GDSET_GCHECKIN_H_ */
