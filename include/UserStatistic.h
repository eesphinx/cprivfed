/*
 * UserStatistic.h
 *
 *  Created on: Jun 23, 2019
 *      
 */

#ifndef PRIV_PLAN_GDSET_USERSTATISTIC_H_
#define PRIV_PLAN_GDSET_USERSTATISTIC_H_

#include <string>

class UserStatistic {
public:

	std::string user_id;
	int count;

	void set(std::string p_user, int p_count);
};




#endif /* PRIV_PLAN_GDSET_USERSTATISTIC_H_ */
