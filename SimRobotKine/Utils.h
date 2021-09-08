/*************************************************
    Copyright:porbots
    Author:dan
    Date:2021-09-08
    Email: liurundan9@163.com
 **************************************************/
 
#ifndef UTILS_H
#define UTILS_H
#include <math.h>
#include <array>

#define PI          3.14159265358979

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


#define DEG2RAD         ( 3.1415926535898/180.0 )
#define RAD2DEG         ( 180.0/3.1415926535898 )

#define IS_ZERO(var)    if(var < 0.0000000001 && var > -0.0000000001){ var = 0;}
#define JUDGE_ZERO(var) ( (var) < 0.0000000001 && (var) > -0.0000000001 )

#define TRF_EPSILON     0.00001

void   print_matrix (double matrix[MATRIX_N][MATRIX_N], int m, int n);
void   print_vector ( std::array<double,6> vec );
double vector_distance ( std::array<double,6> veca, std::array<double,6> vecb );

int matrix_mul(double matrix_a[MATRIX_N][MATRIX_N],
                double matrix_b[MATRIX_N][MATRIX_N],
                double matrix_result[MATRIX_N][MATRIX_N], int m, int n);

int matrix_add(double matrix_a[MATRIX_N][MATRIX_N],
                double matrix_b[MATRIX_N][MATRIX_N],
                double matrix_result[MATRIX_N][MATRIX_N], int m, int n);

void matrix_copy(double matrix_a[MATRIX_N][MATRIX_N],
                 double matrix_b[MATRIX_N][MATRIX_N], int m, int n);

void calculate_matrix_R (double worldrr, double worldrp, double worldry,
                         double (*matrix_R)[MATRIX_N]);

void initmatrix_A       (void* p_param);
void calculate_matrix_A (double matrix[MATRIX_N][MATRIX_N], void *p_param);
void matrix_translate   (double matrix[MATRIX_M][MATRIX_N], int m, int n);

unsigned short DecomposeMatrix( double RM[4][4],
                                double A_actual, double B_actual, double C_actual,
                                double *A, double *B, double *C);
/*
 * ZYX
*/
void MatrixToRPY ( const double RM[4][4], double rpy[3]   );
void RPYToMatrix ( const double rpy[3],   double RM[4][4] );

unsigned short ComposeMatrix(   double RM[4][4],
                                double A, double B, double C);

double Modulo2PI (double actual, double desired);
double ModuloPI (double actual, double desired);
double ValidAngle (double x, double min, double max);
/* declare miscellaneous datatypes and functions */

double RoundToEpsilon(double Value);


#endif // UTILS_H
