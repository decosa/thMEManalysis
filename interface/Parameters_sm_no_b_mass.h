//==========================================================================
// This file has been automatically generated for C++
// MadGraph 5 v. 1.4.8, 2012-07-24
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef Parameters_sm_no_b_mass_H
#define Parameters_sm_no_b_mass_H

#include <complex> 

#include "TopQuarkAnalysis/SingleTop/interface/read_slha.h"
using namespace std; 

class Parameters_sm_no_b_mass
{
  public:

    static Parameters_sm_no_b_mass * getInstance(); 

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double WH, WW, WZ, WT, ymtau, ymt, aS, Gf, aEWM1, MH, MZ, MTA, MT, CKM2x2,
        MZ__exp__2, MZ__exp__4, sqrt__2, MH__exp__2, aEW, MW, sqrt__aEW, ee,
        MW__exp__2, sw2, cw, sqrt__sw2, sw, g1, gw, v, v__exp__2, lam, yt,
        ytau, muH, gw__exp__2, cw__exp__2, ee__exp__2, sw__exp__2;
    std::complex<double> complexi; 
    // Model parameters dependent on aS
    double sqrt__aS, G, G__exp__2; 
    // Model couplings independent of aS
    std::complex<double> GC_17, GC_32, GC_38; 
    // Model couplings dependent on aS


    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 


  private:
    static Parameters_sm_no_b_mass * instance; 
}; 

#endif  // Parameters_sm_no_b_mass_H

