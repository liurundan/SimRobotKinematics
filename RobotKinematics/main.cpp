/*************************************************
    Copyright:porbots
    Author:dan
    Date:2021-09-08
    Email: liurundan9@163.com
 **************************************************/
 
#include <iostream>
#include <stdio.h>
#include "Kinematics.h"

using namespace std;
//6-Axis Robot


double D1 = 399, L1 = 0, L2 = 350, L3 = 42, D3 = 0, D4 = 351, D6 = 82;
/*
 * ABB1200: D1 = 399, L1 = 0, L2 = 350, L3 = 42, D3 = 0, D4 = 351, D6 = 82,
 *
 *                            d       L       a          θ(offset)
*/
//Link_Type DH_Param[] = { { {399, 0,   0}, {-PI/2, 0,     0} },
//                         { {  0, 0, 350}, {    0, 0, -PI/2} },
//                         { {  0, 0,  42}, {-PI/2, 0,     0} },
//                         { {351, 0,   0}, { PI/2, 0,     0} },
//                         { {  0, 0,   0}, {-PI/2, 0,     0} },
//                         { { 82, 0,   0}, {    0, 0,    PI} } };

Link_Type DH_Param[] = { { {  0,  0,    0}, {    0, 0,     0} },
                         { {  0,  0,   D1}, {    0, 0,     0} },
                         { {  0,  0,   L2}, {    0, 0,     0} },
                         { {  D4, 0,   L3}, {    0, 0,     0} },
                         { {  0,  0,    0}, {    0, 0,     0} },
                         { {  D6, 0,    0}, {    0, 0,     0} } };

int main()
{
    cout << "Hello World!" << endl;

    double j[] = { 45, -45, -92, 182, 50, 180 };
    double path[] = { 0, 0, 0, 0, 0, 0};
    double cartpos[] = { 0, 0, 0, 0, 0, 0};

    cout << "input: q =[ ";
    for( int i = 0; i < 6; ++i )
    {
        cout << j[i] << " ";
    }
    cout << "]\n";

    unsigned short ret = ArmDirect(DH_Param,j,path,cartpos);
    cout << "output: cart =[ ";
    for( int i = 0; i < 6; ++i )
    {
        if(i==3) cout << "\n";
        cout << cartpos[i] << " ";
    }
    cout << " ]\n";
    cout << "\n\n";

    double joint[6] ={0};
    ret = ArmInverse( DH_Param, cartpos, j, joint );

    cout << "output: j = [ ";
    for( int i = 0; i < 6; ++i )
    {
        cout << joint[i] << " ";
    }
    cout << "]\n";

    return 0;
}

