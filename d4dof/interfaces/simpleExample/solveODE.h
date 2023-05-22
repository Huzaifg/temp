#ifndef SOLVEODE_H
#define SOLVEODE_H

#include "adept.h"

// Some functions that need to be wrapped
void computeStateJac(std::vector<adept::adouble>& v_st, std::vector<adept::adouble>& old_state, adept::Stack& stack);
void integrate(std::vector<adept::adouble>& v_st, std::vector<adept::adouble>& v_param, std::vector<adept::adouble>& controls, adept::Stack& stack);
void solve(std::vector<adept::adouble>& v_st, std::vector<adept::adouble>& v_params, std::vector<adept::adouble>& controls, std::vector<double> statejac, adept::Stack& stack);


#endif