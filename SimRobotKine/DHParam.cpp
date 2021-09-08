#include "DHParam.h"
#include <string.h>
#include "Utils.h"

DHParam::DHParam()
{
    double D1 = 399, L1 = 0, L2 = 350, L3 = 42, D3 = 0, D4 = 351, D6 = 82;
    offet_j2 = -PI/2;
    offet_j6 = PI;

    /*
     * DH ABB1200-07 robot
     *                 theta,     L,   d,    alpha, negative,   positive
     */
    param_table[0] = { 0,         0,  D1,  -PI/2,    -170,      170      };
    param_table[1] = { 0,        L2,  0,      0,     -100,      135      };
    param_table[2] = { 0,        L3,  0,   -PI/2,    -200,       70      };
    param_table[3] = { 0,         0,  D4,   PI/2,    -270,      270      };
    param_table[4] = { 0,         0,  0,   -PI/2,    -130,      130      };
    param_table[5] = { 0,         0,  D6,      0,    -360,      360      };
}
