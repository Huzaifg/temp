#include <iostream>
#include <fstream>
#include <algorithm>
#include "fourDof.h"
#include <cmath>
#include "adept.h"

using namespace adept;
/// Driver inputs from data file.
/// A driver model based on user inputs provided as time series. If provided as a
/// text file, each line in the file must contain 4 values:
///   time steering throttle braking
/// It is assumed that the time values are unique and soted from 0 to T
void driverInput(std::vector <Entry>& m_data ,const std::string& filename){

    std::ifstream ifile(filename.c_str());
    std::string line;

    
    // get each line
    while(std::getline(ifile,line)){
        std::istringstream iss(line);
        
        double time, steering, throttle, braking;

        // put the stream into our varaibles
        iss >> time >> steering >> throttle >> braking;
        if (iss.fail())
            break;
        // push into our structure
        m_data.push_back(Entry(time, steering, throttle, braking));
        

    }

    ifile.close();
}
/// Function to get the vehicle controls at a given time
/// need to pass the data as well

// template <typename T> // This function is to be used with any type - Stack is optinally null incase we are using the double version of the function
void getControls(std::vector <adouble>& controls, std::vector<Entry>& m_data, const double time, adept::Stack& stack){
   
    
    // if its before time or after time
    if(time <= m_data[0].m_time){
        controls[0] = time;
        controls[1] = m_data[0].m_steering;
        controls[2] = m_data[0].m_throttle;
        controls[3] = m_data[0].m_braking;
        return;
    } else if(time >= m_data.back().m_time){
        controls[0] = time;
        controls[1] = m_data.back().m_steering;
        controls[2] = m_data.back().m_throttle;
        controls[3] = m_data.back().m_braking;
        return;
    }

    // if its within the time, get an iterator and do linear interpolation
    // use compare fucntion earlier defined

    std::vector<Entry>::iterator right =
        std::lower_bound(m_data.begin(), m_data.end(), Entry(time,0,0,0), compareTime); // return first value after Entry
    
    std::vector<Entry>::iterator left = right - 1;

    // linear interplolation
    
    adouble tbar = (time - left->m_time) / (right->m_time - left->m_time);

    controls[0] = time;
    controls[1] = left->m_steering + tbar * (right->m_steering - left->m_steering);
    controls[2] = left->m_throttle + tbar * (right->m_throttle - left->m_throttle);
    controls[3] = left->m_braking + tbar * (right->m_braking - left->m_braking);
}

void getControls(std::vector <double>& controls, std::vector<Entry>& m_data, const double time){
   
    
    // if its before time or after time
    if(time <= m_data[0].m_time){
        controls[0] = time;
        controls[1] = m_data[0].m_steering;
        controls[2] = m_data[0].m_throttle;
        controls[3] = m_data[0].m_braking;
        return;
    } else if(time >= m_data.back().m_time){
        controls[0] = time;
        controls[1] = m_data.back().m_steering;
        controls[2] = m_data.back().m_throttle;
        controls[3] = m_data.back().m_braking;
        return;
    }

    // if its within the time, get an iterator and do linear interpolation
    // use compare fucntion earlier defined

    std::vector<Entry>::iterator right =
        std::lower_bound(m_data.begin(), m_data.end(), Entry(time,0,0,0), compareTime); // return first value after Entry
    
    std::vector<Entry>::iterator left = right - 1;

    // linear interplolation
    
    double tbar = (time - left->m_time) / (right->m_time - left->m_time);

    controls[0] = time;
    controls[1] = left->m_steering + tbar * (right->m_steering - left->m_steering);
    controls[2] = left->m_throttle + tbar * (right->m_throttle - left->m_throttle);
    controls[3] = left->m_braking + tbar * (right->m_braking - left->m_braking);
}

void computeStateJac(VehStates<adouble>& v_st, std::vector<adouble>& old_state, adept::Stack& stack){
    stack.dependent(v_st._x);
    stack.dependent(v_st._y);
    stack.dependent(v_st._theta);
    stack.dependent(v_st._v);
    stack.independent(&old_state[0], 4);
    stack.jacobian(&v_st._stateJac[0]); // Compute and store the State jacobian from the stack
    
    // Once computed, clear dependent and independent as we might need to compute some other jacobian
    stack.clear_dependents();
    stack.clear_independents();

}


void computeInputJac(VehStates<adouble>& v_st, const std::vector<adouble>& controls, adept::Stack& stack){
    stack.dependent(v_st._x);
    stack.dependent(v_st._y);
    stack.dependent(v_st._theta);
    stack.dependent(v_st._v);
    stack.independent(&controls[1],3);
    stack.jacobian(&v_st._inputJac[0]); // Compute and store the input jacobian from the stack

    // Once computed, clear dependent and independent as we might need to compute some other jacobian
    stack.clear_dependents();
    stack.clear_independents();
}

void computeParameterJac(VehStates<adouble>& v_st, const VehParam<adouble>& v_param, adept::Stack& stack){
    stack.dependent(v_st._x);
    stack.dependent(v_st._y);
    stack.dependent(v_st._theta);
    stack.dependent(v_st._v);

    // All parameters that we want sensitivity against
    stack.independent(v_param._l);
    stack.independent(v_param._c0);
    stack.independent(v_param._c1);
    stack.independent(v_param._R); 
    stack.independent(v_param._I);
    stack.independent(v_param._gamma);
    stack.independent(v_param._tau0);
    stack.independent(v_param._omega0);

    stack.jacobian(&v_st._paramJac[0]);

    // Once computed, clear dependent and independent as we might need to compute some other jacobian
    stack.clear_dependents();
    stack.clear_independents();
}

void computeGradYX(const adouble& y, const adouble& x, adouble& out, adept::Stack& stack){
    y.set_gradient(1.0); // Set y as the variable we want the graident for
    stack.compute_adjoint();
    out = x.get_gradient();
}


void computeGradsYXs(const adouble& y, const std::vector<adouble>& xs, std::vector<double>& outs, adept::Stack& stack){
    y.set_gradient(1.0); // Set y as the variable we want the graident for
    stack.compute_adjoint();
    unsigned int lent = xs.size();
    adept::get_gradients(&xs[0], lent, &outs[0]);

}


void integrate(VehStates<adouble>& v_st, const std::vector<adouble>& old, const VehParam<adouble>& v_params, const std::vector <adouble>& controls, adept::Stack& stack){

    v_st._x = old[0] + cos(old[2])*old[3]*v_params._step;
    v_st._y = old[1] + sin(old[2])*old[3]*v_params._step;
    v_st._theta = old[2] + (tan(controls[1]) / v_params._l)*old[3]*v_params._step; // theta update

    adouble rGamma = v_params._R * v_params._gamma;
    adouble f1 = v_params._tau0 * controls[2] - (v_params._tau0 * old[3] / (v_params._omega0 * rGamma));
    adouble f0 = old[3] * v_params._c1 / rGamma + v_params._c0;
    v_st._v = old[3] +  ((rGamma / v_params._I) * (f1 - f0))*v_params._step; // V update
}

// Propogates the state given the controls and the vehicle parameters and also computes and stores the Jacobian information
void solve(VehStates<adouble>& v_st, const VehParam<adouble>& v_params, const std::vector <adouble>& controls, adept::Stack& stack){

    std::vector<adouble> old(4); // vector of current states, independent variables
    old[0] = v_st._x;
    old[1] = v_st._y;
    old[2] = v_st._theta;
    old[3] = v_st._v; // store them
    std::vector<adouble> dparams(5);
    stack.new_recording(); // start a new recording

    integrate(v_st, old, v_params, controls, stack);


    computeStateJac(v_st, old, stack); // Compute and store the Jacobian
    computeInputJac(v_st, controls, stack); // Compute and store the Jacobian with respect to the inputs
    computeParameterJac(v_st, v_params, stack); // Compute and store the Jacobian with respect to the parameters
    computeGradYX(v_st._v, v_params._gamma, v_st._randGrad, stack);
    // Need to add this derivative dependence explicitally -> Also for some reason if I clear gradients here, the next 
    // function gives all zero. Also, when I comment out the above function call, the below one does not work
    double g = 1.;
    v_params._l.add_derivative_dependence(dparams[0], g);
    v_params._c0.add_derivative_dependence(dparams[1], g);
    v_params._c1.add_derivative_dependence(dparams[2], g);
    v_params._R.add_derivative_dependence(dparams[3], g);
    v_params._I.add_derivative_dependence(dparams[4], g);
    computeGradsYXs(v_st._v, dparams, v_st._randGrads, stack);
    // stack.clear_gradients();

}


// Integration for non-AD
void integrate_any(VehStates<double>& v_st, const VehParam<double>& v_params, const std::vector <double>& controls){

    v_st._x = v_st._x + cos(v_st._theta)*v_st._v*v_params._step; // x update
    v_st._y = v_st._y + sin(v_st._theta)*v_st._v*v_params._step; // y update
    v_st._theta = v_st._theta + (tan(controls[1]) / v_params._l) * v_st._v *v_params._step; // theta update

    double rGamma = v_params._R * v_params._gamma;
    double f1 = v_params._tau0 * controls[2] - (v_params._tau0 * v_st._v / (v_params._omega0 * rGamma));
    double f0 = v_st._v * v_params._c1 / rGamma + v_params._c0;
    v_st._v = v_st._v +  ((rGamma / v_params._I) * (f1 - f0))*v_params._step; // V update
}

// Jacobian for non-AD
void jac_any(VehStates<double>& v_st, const VehParam<double>& v_params, const std::vector <double>& controls){

    v_st._stateJac[0] = 1;
    v_st._stateJac[1] = 0;
    v_st._stateJac[2] = 0;
    v_st._stateJac[3] = 0;
    v_st._stateJac[4] = 0;
    v_st._stateJac[5] = 1;
    v_st._stateJac[6] = 0;
    v_st._stateJac[7] = 0;
    v_st._stateJac[8] = -sin(v_st._theta)* v_st._v * v_params._step;
    v_st._stateJac[9] =  cos(v_st._theta) * v_st._v * v_params._step;
    v_st._stateJac[10] = 1;
    v_st._stateJac[11] = 0;
    v_st._stateJac[12] = cos(v_st._theta)*v_params._step;
    v_st._stateJac[13] = sin(v_st._theta)*v_params._step;
    v_st._stateJac[14] = (tan(controls[1]) / v_params._l)*v_params._step;
    v_st._stateJac[15] = 1 - (((v_params._tau0/ v_params._omega0) + v_params._c1) / v_params._I)*v_params._step;
}

// Function for non AD
void solve_any(VehStates<double>& v_st, const VehParam<double>& v_params, const std::vector <double>& controls){
    // compute jacobian anlytically
    jac_any(v_st, v_params, controls);
    integrate_any(v_st, v_params, controls);
}




