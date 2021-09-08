#include "Misc.h"
#include "Frame.h"
#include <stdio.h>



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

void printmatrix(double matrix[MATRIX_N][MATRIX_N], int m, int n)
{
    int i, j;

    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            printf(" %lf ", matrix[i][j]);
       }
    }
    printf(" \n ");
}

void matrix_copy3(double matrix_a[3][3], double matrix_b[3][3])
{
    int i,j;
    for(i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        {
            matrix_b[i][j] = matrix_a[i][j];
        }
    }
}

void printmatrix3( double matrix[3][3] )
{
    printf(" mat = [\n ");
    int i, j;

    for(i=0; i<3; i++)
    {
        for(j=0; j<3; j++)
        {
            printf(" %lf ", matrix[i][j]);
        }
        printf(" \n ");
    }
    printf(" ]\n ");
}
