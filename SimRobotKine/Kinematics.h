/*************************************************
    Copyright:porbots
    Author:dan
    Date:2021-09-08
    Email: liurundan9@163.com
 **************************************************/
#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "DHParam.h"
#include "Utils.h"
#include <vector>

#define TRF_POSE_LEFTY    1
#define TRF_POSE_RIGHTY   0

#define TRF_POSE_BELOW    2 // down eblow
#define TRF_POSE_ABOVE    0 // up eblow

#define TRF_POSE_NONFLIP  4
#define TRF_POSE_FRONT    0

#define TRF_POSE_POSITIVE 0
#define TRF_POSE_NEGATIVE 8

class Kinematics
{
public:
    Kinematics( DHParam* );

    void forward_kine ( const std::array<double,6> joint, std::array<double,6>& cart, char& iFig );
    int  inverse_kine ( const std::array<double,6> cart,  const char iFig, std::array<double,6>& joint );

    double z_offset;

    std::array< std::array<double,6>,8>& get_possible_solution () { return  solution; }
    std::vector< std::array<double, 6>>& get_valid_solution ()    { return  valid_solution; }

private:
    DHParam* dhptr;  
    std::array< std::array<double,6> , 8 >  solution;
    std::vector< std::array<double,6> >     valid_solution;

    double matrix_A[6][MATRIX_N][MATRIX_N];
    double matrix_T[6][MATRIX_N][MATRIX_N];

    void initmatrix_A();
    bool in_range( double val, param_t* par);

    int  get_j456( double jnt[3], double TRPOS[4][4], double sol_1[3], double sol_2[3] );
};

#endif // KINEMATICS_H
