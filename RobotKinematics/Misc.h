/*************************************************
    Copyright:porbots
    Author:dan
    Date:2021-09-08
    Email: liurundan9@163.com
 **************************************************/
 
#ifndef _MISC_TYPES
#define _MISC_TYPES

#include "Frame.h"
#include <string.h>
#include <math.h>

/* declaration of some useful trigonometric functions */
#define TRF_EPSILON     0.00001
#define PI              3.14159265358979

#define sind(x)     sin((x) * PI / 180.0)
#define cosd(x)     cos((x) * PI / 180.0)
#define tand(x)     tan((x) * PI / 180.0)

#define asind(x) 	(asin(x) * 180.0 / PI)
#define acosd(x) 	(acos(x) * 180.0 / PI)
#define atand(x) 	(atan(x) * 180.0 / PI)
#define atan2d(y,x) (atan2(y,x) * 180.0 / PI)

#define sqrt3       sqrt(3.0)
#define sin120      (sqrt(3.0) / 2.0)
#define cos120      -0.5
#define tan60       sqrt(3.0)
#define sin30       0.5
#define tan30       (1.0 / sqrt(3.0))

#define sign(a)     ( (a<0)?-1:1 )

#define MATRIX_1 1
#define MATRIX_M 4
#define MATRIX_N 4

double Modulo2PI(double actual, double desired);
double ModuloPI(double actual, double desired);

/* declare miscellaneous datatypes and functions */

double RoundToEpsilon(double Value);	
void   matrix_copy3 (double matrix_a[3][3], double matrix_b[3][3]);

void   printmatrix (double matrix[MATRIX_N][MATRIX_N], int m, int n);
void   printmatrix3 ( double matrix[3][3] );
#endif
