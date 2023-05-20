#include <iostream>
#include <stdint.h>
#include "fourDof.h"
#include "Timer.h"
#include "StackWrapper.h"

using namespace adept;



int main(int argc, char* argv[]){
    Timer timer;
    timer.print_on_exit(true);
    // input file
    std::string fileName = "acc3.txt";



    std::vector<Entry> driverData;
    driverInput(driverData, fileName);

    
    // First lets time the entire AD code 
    // There can be 2 delays, one is using adoubles and the other is using non-analytical jacobian
    int full_id = timer.new_activity("Simulation with AD");
    // lets initialize a vehicle struct
    // adept::Stack stack;
    StackWrapper stack_wrap;
    VehParam<adouble> veh1_param;
    VehStates<adouble> veh1_st;

    veh1_param._step = 0.001;
    double step = veh1_param._step;
    
    adouble endTime = 10.000;



    // initialize controls to 0
    
    std::vector <adouble> controls(4);

    // now lets run our simulation
    double t = 0;
    int timeStepNo = 0; // time step counter
    timer.start(full_id);
    while(t < (endTime - step/10)){
        
        getControls(controls, driverData, t, stack_wrap.stack);
        solve(veh1_st, veh1_param, controls, stack_wrap.stack);

        t += step;
        timeStepNo += 1;
    }
    timer.stop();
    stack_wrap.stack.pause_recording(); // pause the recording as we are not going to use AD anymore

    std::cout << "With AD" << std::endl;
    std::cout << "states " << std::endl;
    std::cout << veh1_st._x << " ";
    std::cout << veh1_st._y << " ";
    std::cout << veh1_st._theta << " ";
    std::cout << veh1_st._v << " ";
    std::cout << std::endl;

    std::cout<<"Time step "<<timeStepNo<<std::endl;
    std::cout<<"State Jacobian "<<std::endl;
    for(unsigned int i =0; i<16; i++){
        std::cout<<veh1_st._stateJac[i]<<", ";
    }
    std::cout<< std::endl;

    std::cout<<"Input Jacobian "<<std::endl;
    for(unsigned int i =0; i<12; i++){
        std::cout<<veh1_st._inputJac[i]<<", ";
    }
    std::cout<< std::endl;


    std::cout<<"Parameter Jacobian "<<std::endl;
    for(unsigned int i =0; i<32; i++){
        std::cout<<veh1_st._paramJac[i]<<", ";
    }
    std::cout<< std::endl;

    std::cout<<"random grad is "<<veh1_st._randGrad << std::endl;

    std::cout<<"Random Gradients "<<std::endl;
    for(unsigned int i =0; i<10; i++){
        std::cout<<veh1_st._randGrads[i]<<", ";
    }
    std::cout<< std::endl;




    // Now doing the entire simulation without AD
    int no_ad = timer.new_activity("Simulation without AD");

    VehParam<double> vehd_param;
    VehStates<double> vehd_st;

    vehd_param._step = 0.001;
    step = vehd_param._step;

    // initialize controls to 0
    
    std::vector <double> dcontrols(4);

    // now lets run our simulation
    t = 0;
    timeStepNo = 0; // time step counter
    timer.start(no_ad);
    while(t < (endTime - step/10)){
        
        getControls(dcontrols, driverData, t);
        solve_any(vehd_st, vehd_param, dcontrols);

        t += step;
        timeStepNo += 1;
    }
    timer.stop();

    std::cout << "Without AD" << std::endl;
    std::cout << "states " << std::endl;
    std::cout << vehd_st._x << " ";
    std::cout << vehd_st._y << " ";
    std::cout << vehd_st._theta << " ";
    std::cout << vehd_st._v << " ";
    std::cout << std::endl;

    std::cout<<"Time step "<<timeStepNo<<std::endl;
    for(unsigned int i =0; i<16; i++){
        std::cout<<vehd_st._stateJac[i]<<", ";
    }
    std::cout<< std::endl;

    
}

