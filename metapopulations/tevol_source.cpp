/**
 * @file   tevol_source.hpp
 * @brief  ODE system for cloSIR dynamics on a N=2 network
 *
 * Source code.
 * g++ -std=c++11 -O3 -o tevol_source ./tevol_source.cpp $(gsl-config --cflags) $(gsl-config --libs)
 *
 * @author  LHD
 * @since   2020-12-07
 */

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>

#include <boost/multi_array.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "dyn.hpp"

using namespace std;

int main(int argc, const char *argv[]) {
		 
	//Model parameters	
	double lambda = atof(argv[1]);
	double tc = atof(argv[2]);
	double X = atof(argv[3]);
	double Y = atof(argv[4]);
	double rho0 = atof(argv[5]);
	double N1 = atof(argv[6]);
    double epsilon = 1e-3;
    const int dim = 6;
    Sparam param = {lambda, rho0, rho0, X, Y, tc, dim};

    // Integrator parameters
    double t = 0;
    double dt = 1e-6;
    double t_step = 0.2;
    const double eps_abs = 1e-10;
    const double eps_rel = 1e-10;

    // Setting initial conditions
    typedef boost::multi_array<double,2> mat_type;
    typedef mat_type::index index;
    mat_type y(boost::extents[1][dim]);
    fill(y.data(),y.data()+y.num_elements(),0.0);
    // Initial conditions
	y[0][0] = (1.0-epsilon)*N1; //So
	y[0][1] = epsilon*N1; //Io
	y[0][2] = 0.0; //Ro
	y[0][3] = (1.0-epsilon); //Sc
	y[0][4] = epsilon; //Ic
	y[0][5] = 0.0; //Rc

    // Observables
    double lastR = y[0][2]+y[0][5];
    double lastI = y[0][1]+y[0][4];

    // Define GSL odeiv parameters
    const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45;
    gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type, dim);
    gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
    gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (dim);
    gsl_odeiv_system sys = {dydt, NULL, dim, &param};
    
    // Integration
    double diff = 1.0;
    int status(GSL_SUCCESS);
    //for (double t_target = t+t_step; diff > 1e-15; t_target += t_step ) { //stop by difference
    for (double t_target = t+t_step; t_target < 1000; t_target += t_step ) { //stop by time
        while (t < t_target) {
            status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,y.data());
            if (status != GSL_SUCCESS) {
				cout << "SNAFU" << endl;
                break;
			}
        } // end while
        diff = abs(y[0][2]+y[0][5] - lastR);
        lastR = y[0][2]+y[0][5];
        lastI = y[0][1]+y[0][4];
        //cout << t << " " << lambda << " " << tc << " " << rho0 << " " << X << " " << Y << " " << lastR << " " << lastI << "\n";
	} //end while
    //cout.flush();
    cout << lambda << " " << tc << " " << rho0 << " " << X << " " << Y << " " << y[0][2] << " " << y[0][5] << "\n";

    // Free memory
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);
    
    return 0;
}
