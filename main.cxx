#include <iostream>
#include <string.h>
#include <vector>
#include <math.h>
#include <stdexcept>
// #include <python3.8>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;

class MDS{

    private:
    float m; //mass
    float c; //damping coeff
    float k; // spring coeff
    float zeta; //damping ratio
    float wn; //natural frequency
    float time_s[2]; //start time and end time

    vector<float> initial_state; // initial position,velocity;
    bool initial_state_set =false; 


    public:
    
    void print(){cout<<c<<endl;}

    //constructors
    MDS(float m,float c, float k):m(m),c(c),k(k){
        wn = sqrt(k/m);
        zeta = (c/m);
    }

    MDS(float zeta, float wn):zeta(zeta),wn(wn){}

    void set_initial_state(float inx0, float inx0_dot){

        initial_state = {inx0,inx0_dot};
        initial_state_set = true;
    }
    
    vector<float> get_inital_state(){
        return initial_state;
    }

    vector<float> get_config(bool derived = true){

        vector<float> config;
        if(derived){
            config = {wn,zeta};
            }
        else{
            config = {m,c,k};
        }

        return config;
    }

    void set_times(float T, float t0=0){
        time_s[0]=t0;
        time_s[1]=T;
    }
    
    //functions involved in RK4
    float func1(float x, float x_dot){
        return x_dot;
    }

    float func2(float x, float x_dot){
        return (-wn*wn*x - zeta*x_dot);
    }

    //RK4 Solver
    vector<vector<float>> RK4Solver(float T, float t0=0,float dt = 0.001){

        // bool result = false;

        if (initial_state_set = false){
            throw runtime_error("Initial state not state. Use set_initial_state() first.");
        }
        set_times(t0,T);

        vector<float> x_vec; 
        vector<float> xdot_vec;
        vector<float> time_vec;

        float t= t0;
        float x = initial_state[0];
        float x_dot = initial_state[1];

        x_vec.push_back(x);
        xdot_vec.push_back(x_dot);
        time_vec.push_back(t0);

        int ntimes = (T-t0)/dt;

        for(int ii = 0;ii<=ntimes;ii++){
            t= t+dt;

            float k11 = func1(x, x_dot);
            float k12 = func2(x, x_dot);

            float k21 = func1(x + dt / 2 * k11, x_dot + dt / 2 * k12);
            float k22 = func2(x + dt / 2 * k11, x_dot + dt / 2 * k12);

            float k31 = func1(x + dt / 2 * k21, x_dot + dt / 2 * k22);
            float k32 = func2(x + dt / 2 * k21, x_dot + dt / 2 * k22);

            float k41 = func1(x + dt * k31, x_dot + dt * k32);
            float k42 = func2(x + dt * k31, x_dot + dt * k32);

            x = x + dt / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
            x_dot = x_dot + dt / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);

            x_vec.push_back(x);
            xdot_vec.push_back(x_dot);
            time_vec.push_back(t);

        }

    return {x_vec,xdot_vec,time_vec};
    }

};

void plot2D(vector<float> x,vector<float> y, string title = "Mass-Spring-Damper (RK4)", string x_label= "time[s]", string y_label="x-position"){

    // plt::plot(y,x);
    plt::xlabel(x_label);
    plt::ylabel(y_label);
    plt::title("Mass-Spring-Damper (RK4)");
    plt::named_plot("m,c,k= 1,1,0",y,x);
    plt::legend();
    plt::show();

}

int main(){
    MDS myMDS(1,1,0);

    vector<float> config = myMDS.get_config(true);
    // myMDS.print();

    myMDS.set_initial_state(2,0);

    auto solved_RK4 = myMDS.RK4Solver(10);
    plot2D(solved_RK4[0],solved_RK4[2]);

    // auto solved_ana = myMDS.ana_solve(10);
    // plot2D(solved_ana[0],solved_ana[2],"Mass-Spring-Damper (Analytical)");
    

    return 0;
}
