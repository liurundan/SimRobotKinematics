#include "Kinematics.h"

unsigned short ArmDirect(Link_Type Links[6], double JointAxes[6], double PathAxes[6], double Axes[6])
{
    //direct transformations for 6ax robot

        //simplify notation
        double Q1 = JointAxes[0];
        double Q2 = JointAxes[1];
        double Q3 = JointAxes[2];
        double Q4 = JointAxes[3];
        double Q5 = JointAxes[4];
        double Q6 = JointAxes[5];

        double c1 = cosd(Q1);
        double s1 = sind(Q1);
        double c2 = cosd(Q2);
        double s2 = sind(Q2);
        double c3 = cosd(Q3);
        double s3 = sind(Q3);
        double c4 = cosd(Q4);
        double s4 = sind(Q4);
        double c5 = cosd(Q5);
        double s5 = sind(Q5);
        double c6 = cosd(Q6);
        double s6 = sind(Q6);
        double c23 = cosd(Q2+Q3);
        double s23 = sind(Q2+Q3);

        double a1x = Links[1].Offset.X;
        //double a1y = Links[1].Offset.Y;
        double a1z = Links[1].Offset.Z;
        //double a2x = Links[2].Offset.X;
        //double a2y = Links[2].Offset.Y;
        double a2z = Links[2].Offset.Z;
        double a3x = Links[3].Offset.X;
        //double a3y = Links[3].Offset.Y;
        double a3z = Links[3].Offset.Z;
        double a4x = Links[4].Offset.X;
        double a5x = Links[5].Offset.X;

        /* base offset from world frame */
        double ZeroFrame[6];
        ZeroFrame[0] = Links[0].Offset.X;
        ZeroFrame[1] = Links[0].Offset.Y;
        ZeroFrame[2] = Links[0].Offset.Z;
        ZeroFrame[3] = Links[0].Rotation.X;
        ZeroFrame[4] = Links[0].Rotation.Y;
        ZeroFrame[5] = Links[0].Rotation.Z;

        //temporary axes values
        double tmpAxes[6];

        /* check mechanical parameters consistency */
        if ((a1z < 0)||(a2z <= 0)||(a3z < 0)||((a3x+a4x) <= 0)||(a5x <= 0)) {
            return -1;
        }

        double aaa = c3*(a3x+a4x+c5*a5x) + s3*(a3z-c4*s5*a5x);
        double bbb = -s3*(a3x+a4x+c5*a5x) + c3*(a3z-c4*s5*a5x);
        double ddd = a1x + c2*aaa + s2*(bbb+a2z);
        double fff = a1z - s2*aaa + c2*(bbb+a2z);

        /* X axis */
        tmpAxes[0] = c1*ddd - s1*s4*s5*a5x;

        /* Y axis */
        tmpAxes[1] = s1*ddd + c1*s4*s5*a5x;

        /* Z axis */
        tmpAxes[2] = fff;


        /* compose R_total from all joints rotations */
        double R_total[3][3];
        double temp01 = -s1*c4+c1*s23*s4;
        double temp02 = c1*c23*s5+c5*s1*s4+c5*c1*s23*c4;
        double temp11 = c1*c4+s1*s23*s4;
        double temp12 = s1*c23*s5-c1*s4*c5+c5*s1*s23*c4;
        double temp21 = c23*s4;
        double temp22 = -s23*s5+c23*c4*c5;

        R_total[0][0] = c1*c23*c5-s1*s4*s5-c1*s23*c4*s5;
        R_total[0][1] = c6*temp01 + s6*temp02;
        R_total[0][2] = -s6*temp01 + c6*temp02;

        R_total[1][0] = s1*c23*c5+c1*s4*s5-s1*s23*c4*s5;
        R_total[1][1] = c6*temp11 + s6*temp12;
        R_total[1][2] = -s6*temp11 + c6*temp12;

        R_total[2][0] = -s23*c5-c23*c4*s5;
        R_total[2][1] = c6*temp21 + s6*temp22;
        R_total[2][2] = -s6*temp21 + c6*temp22;

        printmatrix3(R_total);
///////
        double rotymat[3][3];//pi/2
        rotymat[0][0] = 0;
        rotymat[0][1] = 0;
        rotymat[0][2] = 1;
        rotymat[1][0] = 0;
        rotymat[1][1] = 1;
        rotymat[1][2] = 0;
        rotymat[2][0] = -1;
        rotymat[2][1] = 0;
        rotymat[2][2] = 0;

        double newdir[3][3]; // set tool Z-axis orentation approach vector. liurundan.
        MatMult(R_total, rotymat, newdir);
//////////////////
        /* A,B,C axis */
        double A,B,C;
        DecomposeMatrix(newdir,PathAxes[3],PathAxes[4],PathAxes[5],&A,&B,&C);
        tmpAxes[3] = A;
        tmpAxes[4] = B;
        tmpAxes[5] = C;

        /* consider zero offset (position and orientation of base point with respect to world origin) */
        SubFrame3D(tmpAxes,ZeroFrame,PathAxes[3],PathAxes[4],PathAxes[5],Axes);

        int i=0;
        for(i=0;i<6;i++)
        {
            Axes[i] = RoundToEpsilon(Axes[i]);
        }

        return 0;
}

unsigned short ArmInverse(Link_Type Links[6], double PathAxes[6], double JointAxes[6], double Axes[6])
{ //inverse transformations for 6ax robot

    /* base offset from world frame */
    double ZeroFrame[6];
    ZeroFrame[0] = Links[0].Offset.X;
    ZeroFrame[1] = Links[0].Offset.Y;
    ZeroFrame[2] = Links[0].Offset.Z;
    ZeroFrame[3] = Links[0].Rotation.X;
    ZeroFrame[4] = Links[0].Rotation.Y;
    ZeroFrame[5] = Links[0].Rotation.Z;

    //temporary axes values
    double tmpAxes[6], WP[6];

    /* consider zero offset (position and orientation of base point with respect to world origin) */
    AddFrame3D(PathAxes,ZeroFrame,PathAxes[3],PathAxes[4],PathAxes[5],tmpAxes);

    //simplify notation
    double X = tmpAxes[0];
    double Y = tmpAxes[1];
    double Z = tmpAxes[2];
    double A = tmpAxes[3];
    double B = tmpAxes[4];
    double C = tmpAxes[5];

    double a1x = Links[1].Offset.X;
    //double a1y = Links[1].Offset.Y;
    double a1z = Links[1].Offset.Z;
    //double a2x = Links[2].Offset.X;
    //double a2y = Links[2].Offset.Y;
    double a2z = Links[2].Offset.Z;
    double a3x = Links[3].Offset.X;
    //double a3y = Links[3].Offset.Y;
    double a3z = Links[3].Offset.Z;
    double a4x = Links[4].Offset.X;
    double a5x = Links[5].Offset.X;

    //keep same pose as current joints configuration
    short Pose = 0;
    if (JointAxes[2] < -90) {
        Pose |= TRF_POSE_CONCAVE; // down elbow
    } else {
        Pose |= TRF_POSE_CONVEX;  // up
    }

    if (JointAxes[4] < 0) {
        Pose |= TRF_POSE_NEGATIVE;
    } else {
        Pose |= TRF_POSE_POSITIVE; // non_flip
    }

    /* check mechanical parameters consistency */
    if ((a1z < 0)||(a2z <= 0)||(a3z < 0)||((a3x+a4x) <= 0)||(a5x <= 0)) {
        return -1;
    }

    /* compose rotation matrix  */
    double R_total[3][3];
    ComposeMatrix(R_total,A,B,C);

    /**************************************************************************
     * @brief rotymat add by liurundan for rotate tool coordinate Z aixs along with tool approach.
     */
    double rotymat[3][3]; //roty(-pi/2) rotate base coordinate Y axis -90.
    rotymat[0][0] = 0;
    rotymat[0][1] = 0;
    rotymat[0][2] = -1;
    rotymat[1][0] = 0;
    rotymat[1][1] = 1;
    rotymat[1][2] = 0;
    rotymat[2][0] = 1;
    rotymat[2][1] = 0;
    rotymat[2][2] = 0;

    double newdir[3][3]; // set tool Z-axis orentation approach vector. liurundan.
    MatMult(R_total, rotymat, newdir);
    matrix_copy3(newdir,R_total);
/*
 * ************************************************************************
*/
    /* calculate wrist point WP from mounting point MP */
    WP[0] = X - a5x * R_total[0][0];
    WP[1] = Y - a5x * R_total[1][0];
    WP[2] = Z - a5x * R_total[2][0];

    /* calculate Q1 */
    /* check for singularity */
    if ((fabs(WP[0]) < TRF_EPSILON)&&(fabs(WP[1]) < TRF_EPSILON))
    {
        Axes[0] = JointAxes[0];
    }
    else
    {
        Axes[0] = atan2d(WP[1],WP[0]);
    }

    /* adjust positions of Q1 with +-PI to bring it closer to desired values */
    Axes[0] = ModuloPI(Axes[0],JointAxes[0]);


    /* consider axis 3 shoulder */
    double b3x = sqrt((a3x+a4x)*(a3x+a4x)+a3z*a3z);

    /* calculate length and height of triangle formed by a2z and b3x */
    double height = WP[2] - a1z;
    double length = sqrt(WP[0]*WP[0] + WP[1]*WP[1]);

    //flip sign of length if Q1 was chosen 180 deg away from atan(Y,X)
    if (fabs(Modulo2PI(Axes[0] - atan2d(WP[1],WP[0]),0))>90)
    {
        length = -length;
    }

    //add a1x correction
    length -= a1x;

    double rho = sqrt(length*length + height*height);

    //check for workspace violations
    if ( (rho > (a2z + b3x + TRF_EPSILON))||(rho < (fabs(a2z-b3x) - TRF_EPSILON)))
    {
        return -1;
    }
    else if (rho > (a2z + b3x)) //adjust impreciseness
    {
        rho = (a2z + b3x);
    }
    else if (rho < fabs(a2z-b3x)) //adjust impreciseness
    {
        rho = fabs(a2z-b3x);
    }

    double alpha = atan2d(height,length);
    double cos_beta = (rho*rho + a2z*a2z - b3x*b3x) / (2.0*a2z*rho);
    double beta = atan2d(sqrt(1-cos_beta*cos_beta),cos_beta);
    double cos_gamma = (a2z*a2z + b3x*b3x - rho*rho) / (2.0*a2z*b3x);
    double gamma = 180.0 - atan2d(sqrt(1-cos_gamma*cos_gamma),cos_gamma);


    if (Pose & TRF_POSE_CONCAVE)
    {
        Axes[1] = 90.0 - alpha + beta;
        Axes[2] = - gamma - atan2d(a3x+a4x,a3z);
    }
    else
    {
        Axes[1] = 90.0 - alpha - beta;
        Axes[2] = gamma - atan2d(a3x+a4x,a3z);
    }

    int i=0;
    for(i=0;i<3;i++)
    {
        Axes[i] = RoundToEpsilon(Axes[i]);
    }

    /* compute wrist rotation matrix */
    /* R_wrist = (R_arm)^T * R_total*/
    double R_arm[3][3];
    double R_wrist[3][3];

    //R_arm = Rz(Q1) * Ry(Q2+Q3)
    double Qy = Axes[1] + Axes[2];
    double Qz = Axes[0];
    R_arm[0][0] = cosd(Qz)*cosd(Qy);
    R_arm[0][1] = -sind(Qz);
    R_arm[0][2] = sind(Qy)*cosd(Qz);
    R_arm[1][0] = cosd(Qy)*sind(Qz);
    R_arm[1][1] = cosd(Qz);
    R_arm[1][2] = sind(Qy)*sind(Qz);
    R_arm[2][0] = -sind(Qy);
    R_arm[2][1] = 0;
    R_arm[2][2] = cosd(Qy);


    //transpose R_arm
    double tmpR;
    tmpR = R_arm[0][1];
    R_arm[0][1] = R_arm[1][0];
    R_arm[1][0] = tmpR;
    tmpR = R_arm[0][2];
    R_arm[0][2] = R_arm[2][0];
    R_arm[2][0] = tmpR;
    tmpR = R_arm[1][2];
    R_arm[1][2] = R_arm[2][1];
    R_arm[2][1] = tmpR;

    MatMult(R_arm,R_total,R_wrist);

    /* extract Q4,Q5,Q6 from wrist rotation matrix as XYX Euler angles */
    /* note that this angle type is not the same as the one used in the decompose matrix function in the frame.h file */

    double A_temp[2],B_temp[2],C_temp[2],ABC_dist[2];
    double A_actual = JointAxes[3];
    double B_actual = JointAxes[4];
    double C_actual = JointAxes[5];

    B_temp[0] = atan2d(sqrt(1-R_wrist[0][0]*R_wrist[0][0]), R_wrist[0][0]);
    B_temp[1] = atan2d(-sqrt(1-R_wrist[0][0]*R_wrist[0][0]), R_wrist[0][0]);

    if (fabs(B_temp[0])>TRF_EPSILON)
    {
        C_temp[0] = atan2d(R_wrist[0][1],R_wrist[0][2]);
        C_temp[1] = atan2d(-R_wrist[0][1],-R_wrist[0][2]);

        A_temp[0] = atan2d(R_wrist[1][0],-R_wrist[2][0]);
        A_temp[1] = atan2d(-R_wrist[1][0],R_wrist[2][0]);
    }
    else
    {	//singularity - choose A=currentQ4
        A_temp[0] = A_temp[1] = A_actual;
        C_temp[0] = C_temp[1] = atan2d(-R_wrist[1][2],R_wrist[2][2]) - A_actual;
    }

    //A, C modulo +-2PI to bring them closer to current values
    A_temp[0] = Modulo2PI(A_temp[0],JointAxes[3]);
    A_temp[1] = Modulo2PI(A_temp[1],JointAxes[3]);

    C_temp[0] = Modulo2PI(C_temp[0],JointAxes[5]);
    C_temp[1] = Modulo2PI(C_temp[1],JointAxes[5]);

    //calculate distance of the two solutions from actual values
    ABC_dist[0] = fabs(A_temp[0]-A_actual) + fabs(B_temp[0]-B_actual) + fabs(C_temp[0]-C_actual);
    ABC_dist[1] = fabs(A_temp[1]-A_actual) + fabs(B_temp[1]-B_actual) + fabs(C_temp[1]-C_actual);

    //keep same pose for wrist
    if (Pose & TRF_POSE_NEGATIVE)
    { //use solution with negative Q5
        if((B_temp[0]<0)&&(B_temp[1]>=0))
        { //use B_temp[0]
            Axes[3] = A_temp[0];
            Axes[4] = B_temp[0];
            Axes[5] = C_temp[0];
        }
        else if((B_temp[1]<0)&&(B_temp[0]>=0))
        { //use B_temp[1]
            Axes[3] = A_temp[1];
            Axes[4] = B_temp[1];
            Axes[5] = C_temp[1];
        }
        else
        { //use closest solution to current values
            if (ABC_dist[0] <= ABC_dist[1])
            {
                Axes[3] = A_temp[0];
                Axes[4] = B_temp[0];
                Axes[5] = C_temp[0];
            }
            else
            {
                Axes[3] = A_temp[1];
                Axes[4] = B_temp[1];
                Axes[5] = C_temp[1];
            }
        }
    }
    else
    { //use solution with positive Q5
        if((B_temp[0]>=0)&&(B_temp[1]<0))
        { //use B_temp[0]
            Axes[3] = A_temp[0];
            Axes[4] = B_temp[0];
            Axes[5] = C_temp[0];
        }
        else if((B_temp[1]>=0)&&(B_temp[0]<0))
        { //use B_temp[1]
            Axes[3] = A_temp[1];
            Axes[4] = B_temp[1];
            Axes[5] = C_temp[1];
        }
        else
        { //use closest solution to current values
            if (ABC_dist[0] <= ABC_dist[1])
            {
                Axes[3] = A_temp[0];
                Axes[4] = B_temp[0];
                Axes[5] = C_temp[0];
            }
            else
            {
                Axes[3] = A_temp[1];
                Axes[4] = B_temp[1];
                Axes[5] = C_temp[1];
            }
        }
    }

    //adjust positions of Q6 with +-2PI to bring it closer to desired value
    Axes[3] = Modulo2PI(Axes[3],JointAxes[3]);
    Axes[5] = Modulo2PI(Axes[5],JointAxes[5]);

    for(i=3;i<6;i++)
    {
        Axes[i] = RoundToEpsilon(Axes[i]);
    }

    return 0;

}

