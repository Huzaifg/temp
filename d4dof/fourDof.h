#ifndef FOURDOF_H
#define FOURDOF_H

#include <stdint.h>
#include <iostream>
#include "adept.h"
/*
Harry's 4 DOF motion model used in state estimation and control. Testing here for sensitivity analysis
*/

// using namespace adept;

// Vehicle Parameters
template <typename T>
struct VehParam{


    // Contructor for the 4 DOF parameters
    VehParam()
        : _l(0.5), _c1(1e-4), _c0(0.02), _R(0.08451952624), _I(1e-3), _gamma(0.334), _tau0(0.09), _omega0(161.185) , _step(1e-3) {}

    
    VehParam(T l, T c1, T c0, T R, T I, T gamma, T tau0, T omega0, double step)
    : _l(l), _c1(c1), _c0(c0), _R(R), _I(I), _gamma(gamma), _tau0(tau0), _omega0(omega0), _step(step) {}


    T _l; // length of car
    T _c1; // Motor resistance that multiplies linearly
    T _c0; // Motor resistance that is constant
    T _R;  // Radius of wheel
    T _I; // Moment of Inertia of wheel
    T _gamma; // Gear Ratio
    T _tau0; // Motor torque at 0 RPM
    T _omega0; // RPM at which motor toruqe goes to 0
    double _step; // Time step used in the integration
};

// Vehicle States
template <typename T>
struct VehStates{
    VehStates()
    : _x(0.), _y(0.), _theta(0.), _v( 0.) {}

    VehStates(T x, T y, T theta, T v)
    : _x(x), _y(y), _theta(theta), _v(v) {}


    T _x; // X position
    T _y; // Y position
    T _theta; // Theta
    T _v; // Speed of the vehicle
    T _randGrad;
    // A array to store the jacobian of the states at each time step
    std::vector<double> _stateJac = std::vector<double>(16);
    std::vector<double> _inputJac = std::vector<double>(12);
    std::vector<double> _paramJac = std::vector<double>(32);
    std::vector<double> _randGrads = std::vector<double>(10);

    


    // Array of adoubles for analytical evaluation
    // std::vector<adouble> _jacAny = std::vector<adouble>(16);
};

// Stores the driver inputs
struct Entry{
    Entry() {} // alias
    // constructor
    Entry(double time, double steering, double throttle, double braking)
        : m_time(time), m_steering(steering), m_throttle(throttle), m_braking(braking) {}

    double m_time;
    double m_steering;
    double m_throttle;
    double m_braking;
};

/// Driver inputs from data file.
void driverInput(std::vector <Entry>& m_data ,const std::string& filename);
/// function needed to compare times for driver input
inline bool compareTime(const Entry& a, const Entry& b){ return a.m_time < b.m_time; };

/// Function to get the vehicle controls at a given time
/// need to pass the data as well
// template <typename T> // This function is to be used with any type - Stack is optinally null incase we are using the double version of the function
void getControls(std::vector <adept::adouble>& controls, std::vector<Entry>& m_data, const double time, adept::Stack& stack);
void getControls(std::vector <double>& controls, std::vector<Entry>& m_data, const double time);
void computeStateJac(VehStates<adept::adouble>& v_st, std::vector<adept::adouble>& old_state, adept::Stack& stack);
void computeInputJac(VehStates<adept::adouble>& v_st, const std::vector<adept::adouble>& controls, adept::Stack& stack);
void computeParameterJac(VehStates<adept::adouble>& v_st, const VehParam<adept::adouble>& v_param, adept::Stack& stack);
void computeGradYX(const adept::adouble& y, const adept::adouble& x, adept::adouble& out, adept::Stack& stack);
void computeGradsYXs(const adept::adouble& y, const std::vector<adept::adouble>& xs, std::vector<double>& outs, adept::Stack& stack);
// This function is to be only used with adouble
void integrate(VehStates<adept::adouble>& v_st, const std::vector<adept::adouble>& old, const VehParam<adept::adouble>& v_params, const std::vector <adept::adouble>& controls, adept::Stack& stack);
void solve(VehStates<adept::adouble>& v_st, const VehParam<adept::adouble>& v_params, const std::vector <adept::adouble>& controls, adept::Stack& stack);

// double version of analytical Jacobian
void integrate_any(VehStates<double>& v_st, const VehParam<double>& v_params, const std::vector <double>& controls);
void jac_any(VehStates<double>& v_st, const VehParam<double>& v_params, const std::vector <double>& controls);
void solve_any(VehStates<double>& v_st, const VehParam<double>& v_params, const std::vector <double>& controls);

#endif
