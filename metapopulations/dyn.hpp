#ifndef DYN_SIMP_HPP_INCLUDED
#define DYN_SIMP_HPP_INCLUDED

#include <boost/multi_array.hpp>

/**
 * @file    dyn_simp.hpp
 * @brief   ODE system for cloSIR dynamics on a N=2 network
 *
 * @author  LHD
 * @since   2020-12-07
 */

struct Sparam {
    const double lambda;
    const double rho012;
    const double rho021;
    const double X;
    const double Y;
    const double tc;
    const int dim;
}; // parameter structure

//********** function dydt definition **************************************************************
int dydt(double t, const double y[], double f[], void * param) {
// ODE system for intra-community growth process

    // Cast parameters
    Sparam& p = *static_cast<Sparam* >(param);

    // Create multi_array reference to y and f
    typedef boost::multi_array_ref<const double,2> CSTmatref_type;
    typedef boost::multi_array_ref<double,2> matref_type;
    typedef CSTmatref_type::index indexref;
    CSTmatref_type yref(y,boost::extents[1][p.dim]);
    matref_type fref(f,boost::extents[1][p.dim]);
    double rho12, rho21, rho22;
    if(t<p.tc) {
        rho12 = p.rho012;
        rho21 = p.rho021;
        rho22 = 1.0;
    } else {
        rho12 = p.rho012*(1.0-p.X);
        rho21 = p.rho021+p.X*p.Y;
        rho22 = (1.0-p.X);
    }

    // Compute derivatives
    //Nodes equations
    //Order: So Io Ro Sc Ic Rc
    fref[0][0] = -p.lambda*yref[0][0]*yref[0][1]-p.lambda*yref[0][0]*(rho12+rho21)*yref[0][4];
    fref[0][1] = +p.lambda*yref[0][0]*yref[0][1]+p.lambda*yref[0][0]*(rho12+rho21)*yref[0][4] - yref[0][1];
    fref[0][2] = +yref[0][1];
    fref[0][3] = -p.lambda*rho22*yref[0][3]*yref[0][4]-p.lambda*yref[0][3]*(rho12+rho21)*yref[0][1];
    fref[0][4] = +p.lambda*rho22*yref[0][3]*yref[0][4]+p.lambda*yref[0][3]*(rho12+rho21)*yref[0][1] - yref[0][4];
    fref[0][5] = +yref[0][4];

    return GSL_SUCCESS;

} //********** end function dydt definition ********************************************************

#endif // DYN_SIMP_HPP_INCLUDED
