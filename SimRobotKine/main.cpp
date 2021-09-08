/*************************************************
	Copyright:porbots
	Author:dan
	Date:2021-09-08
	Email: liurundan9@163.com
 **************************************************/
#include <iostream>
#include <stdio.h>
#include "Kinematics.h"
#include <ctime>
#include <random>
#include <algorithm>
// add by liurundan
using namespace std;

static default_random_engine e(time(0));
static uniform_real_distribution<double> u(-PI, PI);

void generate_random( std::array<double,6>& udata )
{
    for( auto& rd : udata )
    {
        rd = u(e) * RAD2DEG;
    }
}

int main()
{
    cout << "Hello World!" << endl;

    DHParam dh;
    Kinematics kine( &dh );

    uniform_real_distribution<double> u(-PI, PI);
    int test_count = 100;
    int cnt = 0;
    bool inverror = false;
    do
    {
        printf("\n\n*********** count = %d start ***************************\n",cnt+1);

        std::array<double,6> joint = { 30,30,30,0,0,30};
        generate_random(joint);
        joint[4] = 0;

        printf("input_joint:\n");
        print_vector( joint );
        printf("\n");

        std::array<double,6> cart;
        cart.fill(0);
        char iFig = 0;

        printf("------- test forward kinematic ------- \n");

        kine.forward_kine( joint, cart, iFig );
        printf("output_cart:\n");
        print_vector( cart );
        printf(" FIG = 0X%X \n", iFig );
        printf("\n\n");

        std::array<double,6> newjoint;

        printf("------- test inverse kinematic ------- \n");
        kine.inverse_kine( cart, iFig, newjoint );
        std::array< std::array<double,6> , 8 >& sol = kine.get_possible_solution();
        //  std::vector< std::array<double,6> >&sol = kine.get_valid_solution();

        for( unsigned int i = 0; i < sol.size(); ++i )
        {
            printf("sol_%d = [ %9.4f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f]\n",
                   i, sol[i][0], sol[i][1],sol[i][2],sol[i][3],sol[i][4],sol[i][5] );

            std::array<double,6> jntmp = sol[i];
            std::array<double,6> cartmp;
            kine.forward_kine( jntmp, cartmp, iFig );
            double dist = vector_distance( cart, cartmp );
            if( dist > TRF_EPSILON )
            {
                print_vector(cartmp);
                printf("Error inverse kinematic\n");
                inverror = true;
                break;
            }
        }

        if( inverror )
        {
            break;
        }

        if( sol.size() == 0 )
        {
            printf("no valid solution for ikine\n");
        }
        else
            printf("inverse kinematic result success!\n");

        cnt++;
    }
    while( cnt < test_count);
    return 0;
}
