#include "solveODE.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>



using namespace adept;


void computeStateJac(std::vector<adouble>& v_st, std::vector<adouble>& old_state, std::vector<double> stateJac, adept::Stack& stack){
    stack.dependent(v_st[0]);
    stack.dependent(v_st[1]);
    stack.dependent(v_st[2]);
    stack.dependent(v_st[3]);
    stack.independent(&old_state[0], 4);
    stack.jacobian(&stateJac[0]); // Compute and store the State jacobian from the stack
    
    // Once computed, clear dependent and independent as we might need to compute some other jacobian
    stack.clear_dependents();
    stack.clear_independents();

}


// given the olds states, parmaeters and controls, integrate to the next time step using explicit Euler
void integrate(std::vector<adept::adouble>& v_st, std::vector<adept::adouble>& old, std::vector<adept::adouble>& v_param, std::vector<adept::adouble>& controls, adept::Stack& stack){
    v_st[0] = old[0] + cos(old[2])*old[3]*v_param[8];
    v_st[1] = old[1] + sin(old[2])*old[3]*v_param[8];
    v_st[2] = old[2] + (tan(controls[1]) / v_param[0])*old[3]*v_param[8]; // theta update

    adouble rGamma = v_param[3] * v_param[5];
    adouble f1 = v_param[6] * controls[2] - (v_param[6] * old[3] / (v_param[7] * rGamma));
    adouble f0 = old[3] * v_param[1] / rGamma + v_param[2];
    v_st[3] = old[3] +  ((rGamma / v_param[4]) * (f1 - f0))*v_param[8]; // V update
}

// Solve forward and store the state jacobian
void solve(std::vector<adept::adouble>& v_st, std::vector<adept::adouble>& v_params, std::vector<adept::adouble>& controls, std::vector<double> stateJac, adept::Stack& stack){

    std::vector<adouble> old(4);
    old[0] = v_st[0];
    old[1] = v_st[1];
    old[2] = v_st[2];
    old[3] = v_st[3];

    stack.new_recording();

    integrate(v_st, old, v_params, controls, stack);
    computeStateJac(v_st, old, stateJac, stack); // Compute and store the Jacobin
}