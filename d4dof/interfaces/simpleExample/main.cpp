#include <iostream>
#include <stdint.h>
#include "solveODE.h"

using namespace adept;


int main(int argc, char* argv[]){

    Stack stack;

    // Vehicle parameters
    std::vector<adouble> veh1_param(9);
    veh1_param[0] = 0.5;
    veh1_param[1] = 1e-4;
    veh1_param[2] = 0.02;
    veh1_param[3] = 0.08451952624;
    veh1_param[4] = 1e-3;
    veh1_param[5] = 0.334;
    veh1_param[6] = 0.09;
    veh1_param[7] = 161.185;
    veh1_param[8] = 1e-3;


    std::vector<adouble>  veh1_st(4,0.0);
    std::vector<adouble> controls(4,0.0);
    std::vector<double> stateJac(16,0.0);


    adouble t = 0.0;
    while(t < 10.0){

        solve(veh1_st, veh1_param, controls, stateJac, stack);
        t = t + veh1_param[8];
    }



}
