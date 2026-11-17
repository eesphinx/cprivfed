
#include <iostream>
#include <string>
#include <sstream>
#include "GCheckin.h"

// construct a checkin record from a (string) line
GCheckin::GCheckin(const std::string p_data_line) {
    std::istringstream lineStream(p_data_line);
    if (!(lineStream >> user_id >> datetime >> lat >> lng >> loc_id) ) {
    	std::cerr << "unable to read location data from a line" << std::endl;
    	exit(0);
    }
}


//GCheckin::GCheckin() {
//}

//GCheckin::GCheckin(std::string p_user_id, std::string p_datetime, std::string p_lat, std::string p_lng, std::string p_loc_id)
//:user_id(p_user_id), datetime(p_datetime), lat(p_lat), lng(p_lng), loc_id(p_loc_id) {
//}

//void GCheckin::copy_data_from_checkin(const GCheckin& chkin) {
//	user_id = chkin.user_id;
//	datetime = chkin.datetime;
//	lat = chkin.lat;
//	lng = chkin.lng;
//	loc_id = chkin.loc_id;
//}


// construct a checkin record from a (string) line
//void GCheckin::copy_data_from_line(const std::string p_data_line) {
//    std::istringstream lineStream(p_data_line);
//    if (!(lineStream >> user_id >> datetime >> lat >> lng >> loc_id) ) {
//    	std::cerr << "unable to read location data from a line" << std::endl;
//    	exit(0);
//    }
//}
