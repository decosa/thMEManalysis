//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph 5 v. 1.4.8, 2012-07-24
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include <iostream> 
#include <iomanip> 
#include "TopQuarkAnalysis/SingleTop/interface/Parameters_sm_no_b_mass.h"

#define DEBUG  0


// Initialize static instance
Parameters_sm_no_b_mass * Parameters_sm_no_b_mass::instance = 0; 

// Function to get static instance - only one instance per program
Parameters_sm_no_b_mass * Parameters_sm_no_b_mass::getInstance()
{
  if (instance == 0)
    instance = new Parameters_sm_no_b_mass(); 

  return instance; 
}

void Parameters_sm_no_b_mass::setIndependentParameters(SLHAReader& slha)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  WH = slha.get_block_entry("decay", 25, 5.753088e-03); 
  WW = slha.get_block_entry("decay", 24, 2.047600e+00); 
  WZ = slha.get_block_entry("decay", 23, 2.441404e+00); 
  WT = slha.get_block_entry("decay", 6, 1.491500e+00); 
  ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00); 
  ymt = slha.get_block_entry("yukawa", 6, 1.645000e+02); 
  aS = slha.get_block_entry("sminputs", 3, 1.180000e-01); 
  Gf = slha.get_block_entry("sminputs", 2, 1.166390e-05); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.325070e+02); 
  MH = slha.get_block_entry("mass", 25, 1.200000e+02); 
  MZ = slha.get_block_entry("mass", 23, 9.118800e+01); 
  MTA = slha.get_block_entry("mass", 15, 1.777000e+00); 
  MT = slha.get_block_entry("mass", 6, 1.725000e+02); 
  CKM2x2 = 1.; 
  complexi = std::complex<double> (0., 1.); 
  MZ__exp__2 = pow(MZ, 2.); 
  MZ__exp__4 = pow(MZ, 4.); 
  sqrt__2 = sqrt(2.); 
  MH__exp__2 = pow(MH, 2.); 
  aEW = 1./aEWM1; 
  MW = sqrt(MZ__exp__2/2. + sqrt(MZ__exp__4/4. - (aEW * M_PI * MZ__exp__2)/(Gf
      * sqrt__2)));
  sqrt__aEW = sqrt(aEW); 
  ee = 2. * sqrt__aEW * sqrt(M_PI); 
  MW__exp__2 = pow(MW, 2.); 
  sw2 = 1. - MW__exp__2/MZ__exp__2; 
  cw = sqrt(1. - sw2); 
  sqrt__sw2 = sqrt(sw2); 
  sw = sqrt__sw2; 
  g1 = ee/cw; 
  gw = ee/sw; 
  v = (2. * MW * sw)/ee; 
  v__exp__2 = pow(v, 2.); 
  lam = MH__exp__2/(2. * v__exp__2); 
  yt = (ymt * sqrt__2)/v; 
  ytau = (ymtau * sqrt__2)/v; 
  muH = sqrt(lam * v__exp__2); 
  gw__exp__2 = pow(gw, 2.); 
  cw__exp__2 = pow(cw, 2.); 
  ee__exp__2 = pow(ee, 2.); 
  sw__exp__2 = pow(sw, 2.); 
}
void Parameters_sm_no_b_mass::setIndependentCouplings()
{
  GC_17 = (CKM2x2 * ee * complexi)/(sw * sqrt__2); 
  GC_32 = (ee__exp__2 * complexi * v)/(2. * sw__exp__2); 
  GC_38 = -((complexi * yt)/sqrt__2); 
}
void Parameters_sm_no_b_mass::setDependentParameters()
{
  sqrt__aS = sqrt(aS); 
  G = 2. * sqrt__aS * sqrt(M_PI); 
  G__exp__2 = pow(G, 2.); 
}
void Parameters_sm_no_b_mass::setDependentCouplings()
{

}

// Routines for printing out parameters
void Parameters_sm_no_b_mass::printIndependentParameters()
{

  if (DEBUG){
  cout <<  "sm_no_b_mass model parameters independent of event kinematics:" <<
      endl;
  cout << setw(20) <<  "WH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << WH << endl;
  cout << setw(20) <<  "WW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << WW << endl;
  cout << setw(20) <<  "WZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << WZ << endl;
  cout << setw(20) <<  "WT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << WT << endl;
  cout << setw(20) <<  "ymtau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ymtau << endl;
  cout << setw(20) <<  "ymt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ymt << endl;
  cout << setw(20) <<  "aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aS << endl;
  cout << setw(20) <<  "Gf " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << Gf << endl;
  cout << setw(20) <<  "aEWM1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEWM1 << endl;
  cout << setw(20) <<  "MH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MH << endl;
  cout << setw(20) <<  "MZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MZ << endl;
  cout << setw(20) <<  "MTA " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MTA << endl;
  cout << setw(20) <<  "MT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MT << endl;
  cout << setw(20) <<  "CKM2x2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << CKM2x2 << endl;
  cout << setw(20) <<  "complexi " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << complexi << endl;
  cout << setw(20) <<  "MZ__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << MZ__exp__2 << endl;
  cout << setw(20) <<  "MZ__exp__4 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << MZ__exp__4 << endl;
  cout << setw(20) <<  "sqrt__2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sqrt__2 << endl;
  cout << setw(20) <<  "MH__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << MH__exp__2 << endl;
  cout << setw(20) <<  "aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEW << endl;
  cout << setw(20) <<  "MW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << MW << endl;
  cout << setw(20) <<  "sqrt__aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sqrt__aEW << endl;
  cout << setw(20) <<  "ee " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ee << endl;
  cout << setw(20) <<  "MW__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << MW__exp__2 << endl;
  cout << setw(20) <<  "sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sw2 << endl;
  cout << setw(20) <<  "cw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << cw << endl;
  cout << setw(20) <<  "sqrt__sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sqrt__sw2 << endl;
  cout << setw(20) <<  "sw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sw << endl;
  cout << setw(20) <<  "g1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << g1 << endl;
  cout << setw(20) <<  "gw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << gw << endl;
  cout << setw(20) <<  "v " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << v << endl;
  cout << setw(20) <<  "v__exp__2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << v__exp__2 << endl;
  cout << setw(20) <<  "lam " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << lam << endl;
  cout << setw(20) <<  "yt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << yt << endl;
  cout << setw(20) <<  "ytau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ytau << endl;
  cout << setw(20) <<  "muH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << muH << endl;
  cout << setw(20) <<  "gw__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << gw__exp__2 << endl;
  cout << setw(20) <<  "cw__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << cw__exp__2 << endl;
  cout << setw(20) <<  "ee__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << ee__exp__2 << endl;
  cout << setw(20) <<  "sw__exp__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << sw__exp__2 << endl;

  }
}
void Parameters_sm_no_b_mass::printIndependentCouplings()
{

  if (DEBUG){
  cout <<  "sm_no_b_mass model couplings independent of event kinematics:" <<
      endl;
  cout << setw(20) <<  "GC_17 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_17 << endl;
  cout << setw(20) <<  "GC_32 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_32 << endl;
  cout << setw(20) <<  "GC_38 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_38 << endl;
  }
}
void Parameters_sm_no_b_mass::printDependentParameters()
{
  if (DEBUG){
  cout <<  "sm_no_b_mass model parameters dependent on event kinematics:" <<
      endl;
  cout << setw(20) <<  "sqrt__aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << sqrt__aS << endl;
  cout << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G << endl;
  cout << setw(20) <<  "G__exp__2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G__exp__2 << endl;
  }
}
void Parameters_sm_no_b_mass::printDependentCouplings()
{

  if (DEBUG){
  cout <<  "sm_no_b_mass model couplings dependent on event kinematics:" <<
      endl;
  }
}


