//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.4.8, 2012-07-24
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "TopQuarkAnalysis/SingleTop/interface/CPPProcess_bdx_thux.h"
#include "TopQuarkAnalysis/SingleTop/interface/HelAmps_sm_no_b_mass.h"

using namespace MG5_sm_no_b_mass; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: b d~ > t h u~ $$ w+ w- WEIGHTED=6
// Process: b s~ > t h c~ $$ w+ w- WEIGHTED=6

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess_bdx_thux::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm_no_b_mass::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  pars->printIndependentParameters(); 
  pars->printIndependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->MT); 
  mME.push_back(pars->MH); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[1]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void CPPProcess_bdx_thux::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  static bool firsttime = true; 
  if (firsttime)
  {
    pars->printDependentParameters(); 
    pars->printDependentCouplings(); 
    firsttime = false; 
  }

  // Reset color flows
  for(int i = 0; i < 1; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 16; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, 0, -1}, {-1,
      -1, -1, 0, 1}, {-1, -1, 1, 0, -1}, {-1, -1, 1, 0, 1}, {-1, 1, -1, 0, -1},
      {-1, 1, -1, 0, 1}, {-1, 1, 1, 0, -1}, {-1, 1, 1, 0, 1}, {1, -1, -1, 0,
      -1}, {1, -1, -1, 0, 1}, {1, -1, 1, 0, -1}, {1, -1, 1, 0, 1}, {1, 1, -1,
      0, -1}, {1, 1, -1, 0, 1}, {1, 1, 1, 0, -1}, {1, 1, 1, 0, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {36, 36}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_bdx_thux(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_bdx_thux(); 
        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_bdx_thux(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_bdx_thux(); 
      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double CPPProcess_bdx_thux::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 5 && id2 == -1)
  {
    // Add matrix elements for processes with beams (5, -1)
    return matrix_element[0]; 
  }
  else if(id1 == 5 && id2 == -3)
  {
    // Add matrix elements for processes with beams (5, -3)
    return matrix_element[0]; 
  }
  else if(id1 == -1 && id2 == 5)
  {
    // Add matrix elements for processes with beams (-1, 5)
    return matrix_element[1]; 
  }
  else if(id1 == -3 && id2 == 5)
  {
    // Add matrix elements for processes with beams (-3, 5)
    return matrix_element[1]; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void CPPProcess_bdx_thux::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j; 

  // Calculate all wavefunctions
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  sxxxxx(p[perm[3]], +1, w[3]); 
  ixxxxx(p[perm[4]], mME[4], hel[4], -1, w[4]); 
  FFV2_3(w[0], w[2], pars->GC_17, pars->MW, pars->WW, w[5]); 
  FFV2_3(w[4], w[1], pars->GC_17, pars->MW, pars->WW, w[6]); 
  FFS1_1(w[2], w[3], pars->GC_38, pars->MT, pars->WT, w[7]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  VVS1_0(w[6], w[5], w[3], pars->GC_32, amp[0]); 
  FFV2_0(w[0], w[7], w[6], pars->GC_17, amp[1]); 

}
double CPPProcess_bdx_thux::matrix_bdx_thux() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 2; 
  const int ncolor = 1; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1}; 
  static const double cf[ncolor][ncolor] = {{9}}; 

  // Calculate color flows
  jamp[0] = +amp[0] + amp[1]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



