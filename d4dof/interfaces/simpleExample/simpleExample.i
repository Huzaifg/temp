%module simpleExample

%include "std_string.i"
%include "std_vector.i"

%{
#include "adept.h"
#include "solveODE.h"
%}
%include "solveODE.h";

%template(vector_adouble) std::vector <adept::adouble>;
%template(vector_double) std::vector <double>;




