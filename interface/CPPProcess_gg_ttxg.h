//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.4.8, 2012-07-24
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef MG5_Sigma_sm_no_b_mass_gg_ttxg_H
#define MG5_Sigma_sm_no_b_mass_gg_ttxg_H

#include <complex> 
#include <vector> 

#include "TopQuarkAnalysis/SingleTop/interface/Parameters_sm_no_b_mass_ttxj.h"

using namespace std; 

//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > t t~ g WEIGHTED=3
// Process: g g > t t~ g WEIGHTED=3
//--------------------------------------------------------------------------

class CPPProcess_gg_ttxg
{
  public:

    // Constructor.
    CPPProcess_gg_ttxg() {}

    // Initialize process.
    virtual void initProc(string param_card_name); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Info on the subprocess.
    virtual string name() const {return "g g > t t~ g (sm_no_b_mass)";}

    virtual int code() const {return 0;}

    const vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    vector < double * > getMomenta(){return p;}
    void setMomenta(vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 5; 
    static const int nprocesses = 1; 

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 18; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 18; 
    std::complex<double> amp[namplitudes]; 
    double matrix_gg_ttxg(); 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_sm_no_b_mass_ttxj * pars; 

    // vector with external particle masses
    vector<double> mME; 

    // vector with momenta (to be changed each event)
    vector < double * > p; 
    // Initial particle ids
    int id1, id2; 

}; 


#endif  // MG5_Sigma_sm_no_b_mass_gg_ttxg_H
