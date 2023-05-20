%module rom4

%include "std_string.i"
%include "std_vector.i"

%{
#include "fourDof.h"
#include "adept.h"
#include "StackWrapper.h"
%}
%include "fourDof.h";
%include "StackWrapper.h"

%template(vector_adouble) std::vector <adept::adouble>;
%template(vector_double) std::vector <double>;
%template(VehParamD) VehParam<double>;
%template(VehParamAD) VehParam<adept::adouble>;
%template(VehStatesD) VehStates<double>;
%template(VehStatesAD) VehStates<adept::adouble>;
%rename("getControlsAD") getControls(std::vector <adept::adouble>&, std::vector<Entry>&, const double, adept::Stack&);
%rename("getControls") getControls(std::vector <double>&, std::vector<Entry>&, const double, adept::Stack&);



