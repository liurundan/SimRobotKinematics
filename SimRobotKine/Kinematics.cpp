#include "Kinematics.h"

#define   DEBUG     0

Kinematics::Kinematics( DHParam* dhp ):dhptr(dhp)
{
    z_offset = 0;
    for (unsigned int i =0; i < solution.size(); ++i) {
        solution[i].fill(0);
    }
    valid_solution.clear();
}

bool Kinematics::in_range( double val, param_t *par )
{
    if( !par )  return false;
    return (val <= par->maxTheta && val >= par->minTheta );
}

void Kinematics::initmatrix_A()
{
    for ( int i=0; i < 6; ++i)
    {
        calculate_matrix_A( matrix_A[i], &(dhptr->param_table[i]) );
#if DEBUG
        printf("matrix A%d = \n",i+1);
        printmatrix( matrix_A[i], MATRIX_N, MATRIX_N );
#endif
    }
}

void Kinematics::forward_kine( const std::array<double,6> joint, std::array<double,6>& cart, char& iFig)
{
    // prepare input
    for (unsigned int i = 0; i < joint.size(); ++i) {
        dhptr->param_table[i].theta = joint[i] * DEG2RAD;
    }
    dhptr->param_table[1].theta += dhptr->offet_j2;
    dhptr->param_table[5].theta += dhptr->offet_j6;

    //
    initmatrix_A();

    matrix_copy( matrix_A[0], matrix_T[0], MATRIX_N, MATRIX_N );
#if DEBUG
    printf("matrix T%d = \n",0);
    printmatrix( matrix_T[0], MATRIX_N, MATRIX_N );
#endif
    for ( int i = 1; i < 6; ++i )
    {
        matrix_mul( matrix_T[i-1], matrix_A[i], matrix_T[i], MATRIX_N, MATRIX_N );
#if DEBUG
        printf("matrix T%d = \n",i);
        printmatrix( matrix_T[i], MATRIX_N, MATRIX_N );
#endif
    }

    //fetch pos
    for ( int i = 0; i < 3; ++i ) {
        cart[i] = matrix_T[5][i][3];
    }

    //fetch orientation from mat to RPY.
    MatrixToRPY(matrix_T[5], &cart[3]);

    for ( int i = 0; i < 3; ++i ) {
        cart[3+i] *= RAD2DEG;
    }

    iFig = 0;

    double theta1 = dhptr->param_table[0].theta;
    double px_firtjoint = cart[0] * cos( theta1 ) + cart[1] * sin( theta1 );
    bool bArmLefty = ( px_firtjoint >= 0 );
    iFig |= bArmLefty ? TRF_POSE_LEFTY : TRF_POSE_RIGHTY;

    if( bArmLefty )
    {
        if( dhptr->param_table[4].theta > 0 )
        {
            iFig |= TRF_POSE_NONFLIP;
        }

        if( dhptr->param_table[2].theta < 0 )
        {
            iFig |= TRF_POSE_BELOW; // elbow down
        }
    }
    else
    {
        if( dhptr->param_table[4].theta < 0 )
        {
            iFig |= TRF_POSE_NONFLIP;
        }

        if( dhptr->param_table[2].theta > 0 )
        {
            iFig |= TRF_POSE_BELOW; // elbow down
        }
    }
//    printf("pos = [ %.2f, %.2f, %.2f ]\n dir = [ %.2f,%.2f, %.2f ] Fig = %d\n",
//           cart[0],cart[1],cart[2],cart[3],cart[4],cart[5], iFig);

}

int Kinematics::inverse_kine( const std::array<double,6> cart,  const char iFig, std::array<double,6>& joint )
{
    // calculate wirst center position.
    double ex,ey,ez;
    double D6 = dhptr->param_table[5].d;
    double D4 = dhptr->param_table[3].d;
    double D1 = dhptr->param_table[0].d;
    double L1 = dhptr->param_table[0].length;
    double L2 = dhptr->param_table[1].length;
    double L3 = dhptr->param_table[2].length;

    double TR[4][4];
    RPYToMatrix( &cart[3],TR );

    ex = cart[0] - D6 * TR[0][2];
    ey = cart[1] - D6 * TR[1][2];
    ez = cart[2] - D6 * TR[2][2];

    // check if the wrister center is out of range
    double fj = sqrt(ex * ex + ey * ey);
    double pe = fj - L1;
    double qe = ez - D1;
    double rdist = pe*pe + qe*qe;
    double act_dist = sqrt(rdist);
    // max. length from J2 to wrist cen
    double r2 = sqrt(D4 * D4 + L3 * L3);
    double max_dist = r2 + L2;
    if( act_dist > max_dist )
    {
        printf("error: over max workspace \n");
        return -1;
    }

    // calculate theta1
    double theta1,theta2,theta3;
    // J1
    theta1 = atan2( ey, ex );
    double theta1_2 = ValidAngle(theta1+PI, -PI, PI);
    for (unsigned int i = 0; i < solution.size(); ++i) {
        double tmp = i > 3 ? theta1_2 : theta1;
        solution[i][0] = tmp * RAD2DEG;
    }

    // use fig to choose one.
    // calculate theta2
    double r1 = L2;
    double cosaphla = ( r1*r1 + rdist - r2*r2 )/(2*r1*sqrt(rdist));
    double sinaphla = sqrt(1-cosaphla*cosaphla);
    double aphla    = atan2( sinaphla, cosaphla );
    double phi      = atan2(qe,pe);
    theta2          = PI/2 - phi - aphla;
    double theta2_2 = PI/2 - phi + aphla;

    theta2   = ValidAngle( theta2* RAD2DEG, -180, 180 );
    theta2_2 = ValidAngle( theta2_2* RAD2DEG, -180, 180 );

    // calculate theta3
    double cosbeta = ( r1*r1 + r2*r2 - rdist )/(2*r1*r2);
    double sinbeta = sqrt(1-cosbeta*cosbeta);
    double beta    = atan2( sinbeta, cosbeta );
    double angdce  = atan2(D4, L3);
    theta3         = PI - (angdce + beta);
    double theta3_2= -(PI + angdce - beta);

    theta3   = ValidAngle( theta3* RAD2DEG, -180, 180 );
    theta3_2 = ValidAngle( theta3_2* RAD2DEG, -180, 180 );

    //compose solution
    // lefty
    solution[0][1] = solution[1][1] = theta2; // elbow up
    solution[0][2] = solution[1][2] = theta3;

    solution[2][1] = solution[3][1] = theta2_2; // elbow down
    solution[2][2] = solution[3][2] = theta3_2;

    // righty
    solution[4][1] = solution[5][1] = -theta2_2; // elbow up
    solution[4][2] = solution[5][2] = theta3;

    solution[6][1] = solution[7][1] = -theta2; // elbow down
    solution[6][2] = solution[7][2] = theta3_2;

    // calculate theta456
    for ( unsigned int i = 0; i < 4; ++i )
    {
        get_j456( &solution[i*2][0], TR, &solution[i*2][3], &solution[i*2+1][3] );
    }

    valid_solution.clear();

    bool bValid = true;
    for ( unsigned int i = 0; i < solution.size(); ++i)
    {
        bValid = true;
        for (unsigned int j = 0; j < solution[i].size(); ++j )
        {
            if( in_range( solution[i][j], &(dhptr->param_table[j]) ) == false )
            {
                bValid = false;
                break;

            }
        }

        if( bValid )
        {
            valid_solution.push_back(solution[i]);
        }
    }

    return 0;
}

int Kinematics::get_j456( double jnt[3], double TRPOS[4][4], double sol_1[3], double sol_2[3] )
{
    double theta4,theta5,theta6;
    double matrix_a[MATRIX_N][MATRIX_N], matrix_b[MATRIX_N][MATRIX_N];
    double matrix_tmp[MATRIX_N][MATRIX_N];

    dhptr->param_table[0].theta = jnt[0] * DEG2RAD;
    dhptr->param_table[1].theta = (jnt[1]-90)* DEG2RAD;
    dhptr->param_table[2].theta = jnt[2]* DEG2RAD;

/*
 * compute writst rotation matrix
 * Rwritt = (R_arm)^-1 * Rtotal
 * (R_arm)^-1 = (A1*A2*A3)'
 * or use
 * rot(-qy)*rot(-qz)*Rtotal
 */
    calculate_matrix_A(matrix_a, &dhptr->param_table[0]);
    calculate_matrix_A(matrix_b, &dhptr->param_table[1]);
    matrix_mul(matrix_a, matrix_b, matrix_tmp, MATRIX_N, MATRIX_N);

    calculate_matrix_A(matrix_b, &dhptr->param_table[2]);
    matrix_mul(matrix_tmp, matrix_b, matrix_a, MATRIX_N, MATRIX_N);
    matrix_translate(matrix_a, MATRIX_N-1, MATRIX_N-1);

    matrix_mul(matrix_a, TRPOS, matrix_b, MATRIX_N-1, MATRIX_N-1);

    //printmatrix(matrix_b,3,3);

    // calculate j4
    theta4 = atan2( matrix_b[1][2],matrix_b[0][2] );
    double sinq5 = sqrt(matrix_b[1][2]*matrix_b[1][2] + matrix_b[0][2]*matrix_b[0][2]);
    theta5 = atan2(sinq5, matrix_b[2][2]);
    theta6 = atan2( -matrix_b[2][1],matrix_b[2][0] );

    double theta41, theta51, theta61;
    double theta42, theta52, theta62;
    //bool nonFlip = true;
    //if( nonFlip )
    double sind_B = sqrt(1-matrix_b[2][2]*matrix_b[2][2]);
    if( fabs(sind_B) > TRF_EPSILON )
    {
        theta41 = ValidAngle(PI + theta4, -PI, PI);
        theta51 = fabs(theta5);
        theta61 = ValidAngle(PI + theta6, -2*PI, 2*PI);

        theta42 = theta4;
        theta52 = -1* fabs(theta5);
        theta62 = theta6;
    }
    else //singularity q5 = 0
    {
        printf("note: above result lead singularity\n\n");
#if 1
        double actual4 = 0;
        theta41 = actual4;
        theta51 = 0;
        //theta61 = actual4 - sign(-matrix_b[2][0]) * atan2(matrix_b[0][1] * sign(-matrix_b[2][0]),matrix_b[1][1]) + PI;
        theta61 = atan2( matrix_b[1][0], matrix_b[0][0] ) + PI;

        theta42 = -PI; // arbitrary value or actual old value.
        theta52 = 0;
        theta62 = theta61 + PI;
#else

        theta41 = 0;
        theta51 = 0;
        theta61 = atan2( matrix_b[0][1] - matrix_b[1][0], -matrix_b[0][0] - matrix_b[1][1] );
        theta42 = -PI;
        theta52 = 0;
        theta62 = theta61 + PI;
#endif
    }

    sol_1[0] = theta41* RAD2DEG;
    sol_1[1] = theta51* RAD2DEG;
    sol_1[2] = theta61* RAD2DEG;

    sol_2[0] = theta42* RAD2DEG;
    sol_2[1] = theta52* RAD2DEG;
    sol_2[2] = theta62* RAD2DEG;
#if DEBUG
    printf("j456 = %.2f  %.2f %.2f\n %.2f  %.2f %.2f\n",
           theta41* RAD2DEG,theta51* RAD2DEG,theta61* RAD2DEG,
           theta42* RAD2DEG,theta52* RAD2DEG,theta62* RAD2DEG);
#endif
    return 0;
}
