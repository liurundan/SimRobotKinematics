/*************************************************
    Copyright:porbots
    Author:dan
    Date:2021-09-08
    Email: liurundan9@163.com
 **************************************************/
#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "Frame.h"

unsigned short ArmDirect(Link_Type Links[6], double JointAxes[6], double PathAxes[6], double Axes[6]);
unsigned short ArmInverse(Link_Type Links[6], double PathAxes[6], double JointAxes[6], double Axes[6]);
unsigned short ArmWireFrame(Link_Type Links[6], double JointAxes[6], Frame_Type WireFrame[8]);

#endif // KINEMATICS_H
