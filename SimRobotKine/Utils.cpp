#include "Utils.h"
#include <stdio.h>
#include "DHParam.h"
#include <limits>

void calculate_matrix_R(double angle_r, double angle_p, double angle_y,
                        double (*matrix_R)[MATRIX_N])
{
    /*计算旋转矩阵 */
    double r_c, r_s, p_c, p_s, y_c, y_s;

    angle_r *= DEG2RAD;
    angle_p *= DEG2RAD;
    angle_y *= DEG2RAD;

    r_c = cos( angle_r );
    IS_ZERO(r_c);
    r_s = sin( angle_r );
    IS_ZERO(r_s);
    p_c = cos( angle_p );
    IS_ZERO(p_c);
    p_s = sin( angle_p );
    IS_ZERO(p_s);
    y_c = cos( angle_y );
    IS_ZERO(p_c);
    y_s = sin( angle_y );
    IS_ZERO(y_s);

    matrix_R[0][0] = r_c * p_c;
    matrix_R[0][1] = r_c * p_s * y_s - r_s * y_c;
    matrix_R[0][2] = r_c * p_s * y_c + r_s * y_s;

    matrix_R[1][0] = r_s * p_c;
    matrix_R[1][1] = r_s * p_s * y_s + r_c * y_c;
    matrix_R[1][2] = r_s * p_s * y_c - r_c * y_s;

    matrix_R[2][0] = -p_s;
    matrix_R[2][1] = p_c * y_s;
    matrix_R[2][2] = p_c * y_c;
}

//根据关节参数计算T矩阵
void calculate_matrix_A(double matrix[MATRIX_N][MATRIX_N], void *p_param)
{
    param_t *dh = (param_t*)p_param;
    double *pmatrix=(double *)matrix;
    double  var_c, var_s, angle_c, angle_s;

    var_c = cos(dh->theta);
    IS_ZERO(var_c);
    var_s = sin(dh->theta);
    IS_ZERO(var_s);
    angle_c = cos(dh->alpha);
    IS_ZERO(angle_c);
    angle_s = sin(dh->alpha);
    IS_ZERO(angle_s);
    *pmatrix++ = var_c;
    *pmatrix++ = -var_s * angle_c;
    *pmatrix++ = var_s * angle_s;
    *pmatrix++ = dh->length * var_c;

    *pmatrix++ = var_s;
    *pmatrix++ = var_c * angle_c;
    *pmatrix++ = -var_c *angle_s;
    *pmatrix++ = dh->length * var_s;

    *pmatrix++ =0;
    *pmatrix++ = angle_s;
    *pmatrix++ = angle_c;
    *pmatrix++ = dh->d;

    *pmatrix++ =0;
    *pmatrix++ =0;
    *pmatrix++ =0;
    *pmatrix =1;
}

void matrix_copy(double matrix_a[MATRIX_N][MATRIX_N],
                 double matrix_b[MATRIX_N][MATRIX_N], int m, int n)
{
    int i,j;
    for(i=0; i<m; i++)
    {
        for(j=0; j<n; j++)
        {
            matrix_b[i][j] = matrix_a[i][j];
        }
    }
}

int matrix_mul(double matrix_a[MATRIX_N][MATRIX_N],
               double matrix_b[MATRIX_N][MATRIX_N],
               double matrix_result[MATRIX_N][MATRIX_N], int m, int n)
{
    int i,j,k;
    double sum;
    double matrix_tmp[MATRIX_N][MATRIX_N]={{0},{0},{0},{0}};

    /*嵌套循环计算结果矩阵的每个元素*/
    for(i=0; i<m; i++)
    {
        for(j=0; j<n; j++)
        {
            /*按照矩阵乘法的规则计算结果矩阵的i*j元素*/
            sum=0;
            for(k=0; k<n; k++)
                sum += matrix_a[i][k] * matrix_b[k][j];
            matrix_tmp[i][j] = sum;
        }
    }
    matrix_copy(matrix_tmp, matrix_result, MATRIX_N, MATRIX_N);
    return 0;
}

int matrix_add(double matrix_a[MATRIX_N][MATRIX_N],
            double matrix_b[MATRIX_N][MATRIX_N],
            double matrix_result[MATRIX_N][MATRIX_N], int m, int n)
{
    int i,j;
    double matrix_tmp[MATRIX_N][MATRIX_N]={{0},{0},{0},{0}};

    for(i=0; i<m; i++)
    {
        for(j=0; j<n; j++)
        {
            matrix_tmp[i][j] = matrix_a[i][j] + matrix_b[i][j];
        }
    }
    matrix_copy(matrix_tmp, matrix_result, MATRIX_N, MATRIX_N);
    return 0;
}

//矩阵转置
void matrix_translate(double matrix[MATRIX_M][MATRIX_N], int m, int n)
{
    double m_tmp;
    int i, j, k;

    for(i=0, j=0; i<m; i++, j++)
    {
        for(k=j; k<n; k++)
        {
            if(i == k) continue;
            m_tmp = matrix[i][k];
            matrix[i][k] = matrix[k][i];
            matrix[k][i] = m_tmp;
        }
    }
}

void print_matrix(double matrix[MATRIX_N][MATRIX_N], int m, int n)
{
    int i, j;

    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            printf(" %lf ", matrix[i][j]);
       }
       printf(" \n ");
    }
    printf(" \n ");
}

unsigned short ComposeMatrix(double RM[4][4], double A, double B, double C)
{// compose rotation matrix from RPY angles

    RM[0][0] = cosd(C)*cosd(B);
    RM[0][1] = cosd(C)*sind(B)*sind(A) - sind(C)*cosd(A);
    RM[0][2] = cosd(C)*sind(B)*cosd(A) + sind(C)*sind(A);
    RM[1][0] = sind(C)*cosd(B);
    RM[1][1] = sind(C)*sind(B)*sind(A) + cosd(C)*cosd(A);
    RM[1][2] = sind(C)*sind(B)*cosd(A) - cosd(C)*sind(A);
    RM[2][0] = - sind(B);
    RM[2][1] = cosd(B)*sind(A);
    RM[2][2] = cosd(B)*cosd(A);

    return 0;
}

unsigned short DecomposeMatrix(double RM[4][4], double A_actual, double B_actual, double C_actual, double *A, double *B, double *C)
{// decompose rotation matrix into RPY angles, ZYX

    double A_temp[2],B_temp[2],C_temp[2],ABC_dist[2];
    double cosd_B = sqrt(1-RM[2][0]*RM[2][0]);

    B_temp[0] = atan2d(-RM[2][0],cosd_B);
    B_temp[1] = atan2d(-RM[2][0],-cosd_B);

    if (fabs(cosd_B)>TRF_EPSILON) {
        C_temp[0] = atan2d(RM[1][0],RM[0][0]);
        C_temp[1] = atan2d(-RM[1][0],-RM[0][0]);

        A_temp[0] = atan2d(RM[2][1],RM[2][2]);
        A_temp[1] = atan2d(-RM[2][1],-RM[2][2]);
    } else {	//singularity - choose A=A_actual
        A_temp[0] = A_temp[1] = A_actual;
        C_temp[0] = C_temp[1] = A_actual - sign(-RM[2][0]) * atan2d(RM[0][1]*sign(-RM[2][0]),RM[1][1]);
    }

    //A, C modulo +-2PI to bring them closer to current values
    A_temp[0] = Modulo2PI(A_temp[0],A_actual);
    A_temp[1] = Modulo2PI(A_temp[1],A_actual);

    C_temp[0] = Modulo2PI(C_temp[0],C_actual);
    C_temp[1] = Modulo2PI(C_temp[1],C_actual);

    //calculate distance of the two solutions from actual values
    ABC_dist[0] = fabs(A_temp[0]-A_actual) + fabs(B_temp[0]-B_actual) + fabs(C_temp[0]-C_actual);
    ABC_dist[1] = fabs(A_temp[1]-A_actual) + fabs(B_temp[1]-B_actual) + fabs(C_temp[1]-C_actual);

    //keep same pose for wrist
    if (B_actual < 0) {//use solution with negative B
        if((B_temp[0]<0)&&(B_temp[1]>=0)) { //use B_temp[0]
            *A=A_temp[0];
            *B=B_temp[0];
            *C=C_temp[0];
        } else if((B_temp[1]<0)&&(B_temp[0]>=0)) { //use B_temp[1]
            *A=A_temp[1];
            *B=B_temp[1];
            *C=C_temp[1];
        } else { //use closest solution to current values
            if (ABC_dist[0] <= ABC_dist[1]) {
                *A = A_temp[0];
                *B = B_temp[0];
                *C = C_temp[0];
            } else {
                *A = A_temp[1];
                *B = B_temp[1];
                *C = C_temp[1];
            }
        }
    } else {//use solution with positive B
        if((B_temp[0]>=0)&&(B_temp[1]<0)) { //use B_temp[0]
            *A = A_temp[0];
            *B = B_temp[0];
            *C = C_temp[0];
        } else if((B_temp[1]>=0)&&(B_temp[0]<0)) { //use B_temp[1]
            *A = A_temp[1];
            *B = B_temp[1];
            *C = C_temp[1];
        } else { //use closest solution to current values
            if (ABC_dist[0] <= ABC_dist[1]) {
                *A = A_temp[0];
                *B = B_temp[0];
                *C = C_temp[0];
            } else {
                *A = A_temp[1];
                *B = B_temp[1];
                *C = C_temp[1];
            }
        }
    }

    //adjust positions of C with +-2PI to bring it closer to desired value
    *C = Modulo2PI(*C,C_actual);

    return 0;
}

/*
 *
 *
*/

double RoundToEpsilon(double Value)
{
    int nValue = round(Value);
    if (fabs(Value - nValue) < TRF_EPSILON)
        return nValue;
    else
        return Value;
}

double Modulo2PI(double actual, double desired)
{// add or subtract multiples of 2PI to the actual angle to bring it closer to the desired value
    double delta = (actual - desired);
    int n = ((int)delta + 180 * sign(delta)) / 360;
    return (actual - (360.0 * n));
}

double ModuloPI(double actual, double desired)
{// add or subtract multiples of PI to the actual angle to bring it closer to the desired value
    actual = Modulo2PI(actual,desired);
    double delta = (actual - desired);
    if (delta > 90)		actual -=180;
    if (delta < -90)	actual +=180;
    return actual;
}

// ZYX
void MatrixToRPY( const double RM[4][4], double rpy[3])
{
    if (RM[2][0]==1) {
         rpy[0] = atan2(-RM[0][1],-RM[0][2]);
         rpy[1] = -M_PI/2;
         rpy[2] = 0.0;
      } else if (RM[2][0]==-1) {
        rpy[0] = atan2(RM[0][1],RM[0][2]);
        rpy[1] = M_PI/2;
        rpy[2] = 0.0;
     }else {
        rpy[0] = atan2(RM[2][1], RM[2][2]);
        rpy[1] = atan2(-RM[2][0], sqrt(RM[0][0]*RM[0][0] + RM[1][0]*RM[1][0]));
        rpy[2] = atan2(RM[1][0], RM[0][0]);
     }
}

void RPYToMatrix(const double rpy[3], double RM[4][4])
{
    double A = rpy[0];
    double B = rpy[1];
    double C = rpy[2];

    RM[0][0] = cosd(C)*cosd(B);
    RM[0][1] = cosd(C)*sind(B)*sind(A) - sind(C)*cosd(A);
    RM[0][2] = cosd(C)*sind(B)*cosd(A) + sind(C)*sind(A);
    RM[1][0] = sind(C)*cosd(B);
    RM[1][1] = sind(C)*sind(B)*sind(A) + cosd(C)*cosd(A);
    RM[1][2] = sind(C)*sind(B)*cosd(A) - cosd(C)*sind(A);
    RM[2][0] = -sind(B);
    RM[2][1] = cosd(B)*sind(A);
    RM[2][2] = cosd(B)*cosd(A);
}

double ValidAngle(double x, double min, double max)
{
    double  modulus = max - min;
    x -= (modulus * floor(x/modulus));
    x = fmod(x - min, modulus) + min;
    return ( fabs(x - min) >=  TRF_EPSILON ) ? x: max;
}

void print_vector(std::array<double, 6> vec)
{
    printf(" \n[ ");
    for(unsigned int i=0; i<vec.size(); i++)
    {
        printf(" %lf ", vec[i]);
    }
    printf(" ]\n ");
}

double vector_distance(std::array<double, 6> veca, std::array<double, 6> vecb)
{
    double dist = 0;
    if( veca.size() != vecb.size() )
        return -1;

    for(unsigned int i = 0; i < veca.size(); ++i )
    {
        dist += pow( (veca[i]-vecb[i]), 2 );
    }
    return dist;
}
