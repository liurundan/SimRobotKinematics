/*************************************************
    Copyright:porbots
    Author:dan
    Date:2021-09-08
    Email: liurundan9@163.com
 **************************************************/
 
#ifndef FRAME_H
#define FRAME_H

#include <math.h>
#include "Misc.h"

#define TRF_POSE_BACK     4
#define TRF_POSE_FRONT    0

#define TRF_POSE_CONCAVE  2 // down eblow
#define TRF_POSE_CONVEX   0 // up eblow

#define TRF_POSE_LEFT     1
#define TRF_POSE_RIGHT    0

#define TRF_POSE_POSITIVE 0
#define TRF_POSE_NEGATIVE 8

#define TRF_DIRECT        0
#define TRF_INVERSE       1

typedef struct Coord_Type
{	double X;
    double Y;
    double Z;
} Coord_Type;

typedef struct Link_Type
{	struct Coord_Type Offset;
    struct Coord_Type Rotation;
} Link_Type;

typedef struct Frame_Type
{	double Axes[6];
} Frame_Type;

/* quaternion datatype */
typedef struct Quat_Type
{	
	double w,x,y,z;
} Quat_Type;

/* Declaration of common functions used for frames operations: translations, rotations... */

unsigned short ComposeMatrix(double RM[3][3], double A, double B, double C);
unsigned short DecomposeMatrix(double RM[3][3], double A_actual, double B_actual, double C_actual, double *A, double *B, double *C);
unsigned short MatMult(double M1[3][3],double M2[3][3],double M3[3][3]);
unsigned short SubFrame3D(double P1[6],double F1[6], double A_actual, double B_actual, double C_actual, double P0[6]);
unsigned short AddFrame3D(double P1[6],double F1[6], double A_actual, double B_actual, double C_actual, double P0[6]);
unsigned short SubFrame2D(double P1[6],double F1[6], double P0[6]);
unsigned short AddFrame2D(double P1[6],double F1[6], double P0[6]);
unsigned short SubFrameTool2D(double Path[6], double Frame[6], double Tool[6], double Mount[6]);
unsigned short AddFrameTool2D(double Mount[6], double Frame[6], double Tool[6], double Path[6]);
unsigned short AddFrameTool3D(double Mount[6], double Frame[6], double Tool[6], double A_actual, double B_actual, double C_actual, double Path[6]);
unsigned short SubFrameTool3D(double Path[6], double Frame[6], double Tool[6], double A_actual, double B_actual, double C_actual, double Mount[6]);
unsigned short NormalizeQuat(Quat_Type* q);
unsigned short MatrixToQuat(double RM[3][3], Quat_Type* q);
unsigned short QuatToMatrix(Quat_Type q, double RM[3][3]);
double AngleBetweenQuat(Quat_Type q1, Quat_Type* q2);
unsigned short Slerp(Quat_Type q1, Quat_Type q2, Quat_Type* q, double angle, double u);
unsigned short EulerToQuat(double A, double B, double C, Quat_Type* q);
unsigned short QuatToEuler(Quat_Type q, double A_actual, double B_actual, double C_actual, double *A, double *B, double *C);

#endif

