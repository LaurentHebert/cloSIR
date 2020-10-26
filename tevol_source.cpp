/**
 * @file   tevol_source.hpp
 * @brief  ODE system for gatherings following SIR dynamics with grathering closure (cloSIR model)
 *
 * Source code.
 * Compilation command :
 * g++ -std=c++11 -O3 -o tevol_source ./tevol_source.cpp $(gsl-config --cflags) $(gsl-config --libs)
 * Example call :
 * ./tevol_source.cpp 1.5 10 0.75 0.5
 *
 * @author  LHD
 * @since   2020-03-01
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
    double epsilon = 1e-3;
    const int dim = 6;
    Sparam param = {lambda,dim};

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
	y[0][0] = (1.0-epsilon); //So
	y[0][1] = epsilon; //Io
	y[0][2] = 0.0; //Ro
	y[0][3] = 0.0; //Sc
	y[0][4] = 0.0; //Ic
	y[0][5] = 0.0; //Rc

    // Define GSL odeiv parameters
    const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45;
    gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type, dim);
    gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
    gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (dim);
    gsl_odeiv_system sys = {dydt, NULL, dim, &param};
	
	//Integration up to intervention
    int status(GSL_SUCCESS);
    for (double t_target = t+t_step; t_target < tc; t_target += t_step ) { 
        while (t < t_target) {
            status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,y.data());
            if (status != GSL_SUCCESS) {
				cout << "SNAFU" << endl;
                break;
			}
        }
        double currentR = y[0][2];
        double currentI = y[0][1];
        cout << t << ", " << lambda << ", " << tc << ", " << X << ", " << Y << ", " << currentR << ", " << currentI << "\n";
    }

    // Redistribute churchgoers
    double Sgoers = y[0][0];
    double Igoers = y[0][1];
    double Rgoers = y[0][2];
    y[0][0] = Sgoers*(1.0+X*Y/(1.0-X)); //So
	y[0][1] = Igoers*(1.0+X*Y/(1.0-X)); //Io
	y[0][2] = Rgoers*(1.0+X*Y/(1.0-X)); //Ro
	y[0][3] = Sgoers*(1.0-Y); //Sc
	y[0][4] = Igoers*(1.0-Y); //Ic
	y[0][5] = Rgoers*(1.0-Y); //Rc
    double lastR = (1.0-X)*y[0][2]+X*y[0][5];
    double lastI = (1.0-X)*y[0][1]+X*y[0][4];

    
    double diff = 1.0;
    //for (double t_target = t+t_step; diff > 1e-15; t_target += t_step ) { //stop by difference
    for (double t_target = t+t_step; t_target < 200; t_target += t_step ) { //stop by time
        while (t < t_target) {
            status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,y.data());
            if (status != GSL_SUCCESS) {
				cout << "SNAFU" << endl;
                break;
			}
        } // end while
        diff = abs((1.0-X)*y[0][2]+X*y[0][5] - lastR);
        lastR = (1.0-X)*y[0][2]+X*y[0][5];
        lastI = (1.0-X)*y[0][1]+X*y[0][4];
        cout << t << ", " << lambda << ", " << tc << ", " << X << ", " << Y << ", " << lastR << ", " << lastI << "\n";
	} //end while
    //cout << lambda << " " << tc << " " << X << " " << Y << " " << lastR << " " << y[0][2] << " " << y[0][5] << "\n";

    cout.flush();

    // Free memory
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);
    
    return 0;
}
