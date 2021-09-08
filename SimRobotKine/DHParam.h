/*************************************************
    Copyright:porbots
    Author:dan
    Date:2021-09-08
    Email: liurundan9@163.com
 **************************************************/
 
#ifndef DHPARAM_H
#define DHPARAM_H
#include <array>

typedef struct {
    double theta;  //joint variable
    double length;
    double d;
    double alpha;

    double minTheta; // soft limit negative, in degree
    double maxTheta;
}param_t;

class DHParam
{
public:
    DHParam();

    std::array<param_t,6> param_table;

    double offet_j2;
    double offet_j6;

};

#endif // DHPARAM_H
