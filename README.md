# cloSIR

Implementation of the cloSIR model Althouse et al. (2020)

### Codes

Integration of the ODEs is done in C++ using the Gnu Scientific Library and Boost library. 
Compilation line and sample call are in the header of tevol_source.cpp.

The code tevol_source.cpp takes in 4 parameters: 
1. normalized transmission rate (lambda)
2. time of intervention (t_c)
3. scale of intervention (X), 
4. and fraction of non-compliance (Y) 
