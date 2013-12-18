#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>
#include <algorithm>

#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TMatrixT.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/AllIntegrationTypes.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "TopQuarkAnalysis/SingleTop/interface/MEIntegratorNew.h"
#include "TopQuarkAnalysis/SingleTop/interface/Samples.h"
#include "TopQuarkAnalysis/SingleTop/interface/MECombination.h"

 




#define GENJETDR  0.3
#define MAX_REEVAL_TRIES 3
#define PI TMath::Pi()
#define NMAXPERMUT 30
#define NMAXMASS   20
#define NMAXJETS   10

using namespace std;


typedef struct 
{
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;


typedef struct 
{
  int run;
  int lumi;
  int event;
  int json;
} EventInfo;
 

//struct sorterByPt {
//bool operator() (float i,float j) const { return (i>j);}
//};


typedef struct
{
  TLorentzVector p4;
  float csv;
} JetObservable;

struct JetObservableListerByPt
{
  bool operator()( const JetObservable &lx, const JetObservable& rx ) const {
    return (lx.p4).Pt() > (rx.p4).Pt();
  }
};

struct JetObservableListerByCSV
{
  bool operator()( const JetObservable &lx, const JetObservable& rx ) const {
    return (lx.csv) > (rx.csv);
  }
};


typedef struct
{
  float bmass;
  float bpt;
  float beta;
  float bphi;
  float bstatus;
  float wdau1mass;
  float wdau1pt;
  float wdau1eta;
  float wdau1phi;
  float wdau1id;
  float wdau2mass;
  float wdau2pt;
  float wdau2eta;
  float wdau2phi;
  float wdau2id;
} genTopInfo;

typedef struct 
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
} genParticleInfo;

typedef struct 
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid1;
  float momid2;
} genParton;


//https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
float resolutionBias(float eta, int shift)
{
  if(eta < 0.5){
    // bias 1s up
    if(shift ==+1) return 0.058;
    // nominal bias
    if(shift == 0) return 0.052;  
    // bias 1s down
    if(shift ==-1) return 0.048;
  }
  if(eta < 1.1){
    // bias 1s up
    if(shift ==+1) return 0.063;
    // nominal bias
    if(shift == 0) return 0.057;    
    // bias 1s down
    if(shift ==-1) return 0.051;
  }
  if(eta < 1.7){
    // bias 1s up
    if(shift ==+1) return 0.102;
    // nominal bias
    if(shift == 0) return 0.096;    
    // bias 1s down
    if(shift ==-1) return 0.090;
  }
  if(eta < 2.3){
    // bias 1s up
    if(shift ==+1) return 0.224;


    // nominal bias
    if(shift == 0) return 0.134;    
    // bias 1s down
    if(shift ==-1) return 0.044;
  }
  if(eta < 5.0){
    // bias 1s up
    if(shift ==+1) return 0.488;
    // nominal bias
    if(shift == 0) return 0.288;    
    // bias 1s down
    if(shift ==-1) return 0.088;
  }
  return 0;
}


int main(int argc,  char* argv[])
{

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@ FWLITE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

 
  //std::cout << "MEAnalysisNew" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();


  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@ CONFIGURATION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */



  //  PsetReader reader("process");
  //  in = reader.read(argc, argv);


  PythonProcessDesc builder(argv[1],argc, argv);
    const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  // SAMPLES
  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;

  std::vector< std::string > inFileNames( in.getParameter< std::vector<std::string> >  ("inFileNames" ) );

  std::string outFileName( in.getParameter<std::string>  ("outFileName" ) );
  std::string pathToFile ( in.getParameter<std::string>  ("pathToFile" ) );
  std::string ordering   ( in.getParameter<std::string>  ("ordering" ) );
  std::string pathToTF   ( in.getParameter<std::string>  ("pathToTF") );
  std::string pathToCP   ( in.getParameter<std::string>  ("pathToCP") );
  bool   verbose         ( in.getParameter<bool>         ("verbose" ) );
  
  // PARAMETERS
  double lumi               ( in.getParameter<double>("lumi") );
  float  MH                 ( in.getUntrackedParameter<double> ("MH",     125.));
  float  MT                 ( in.getUntrackedParameter<double> ("MT",    174.3));
  float  MW                 ( in.getUntrackedParameter<double> ("MW",    80.19));
  float  MwL                ( in.getUntrackedParameter<double> ("MwL",      60));
  float  MwH                ( in.getUntrackedParameter<double> ("MwH",     100));
  float  btag_prob_cut_6jets( in.getUntrackedParameter<double> ("btag_prob_cut_6jets",    0.));
  float  btag_prob_cut_5jets( in.getUntrackedParameter<double> ("btag_prob_cut_5jets",    0.));
  float  btag_prob_cut_4jets( in.getUntrackedParameter<double> ("btag_prob_cut_4jets",    0.));
  double maxChi2_           ( in.getUntrackedParameter<double> ("maxChi2", 2.5));
  vector<double> massesH    ( in.getParameter<vector<double> > ("massesH"));
  vector<double> massesT    ( in.getParameter<vector<double> > ("massesT"));

  // FLAGS
  int   switchoffOL      ( in.getUntrackedParameter<int>    ("switchoffOL",     0));
  int   speedup          ( in.getUntrackedParameter<int>    ("speedup",         0));
  int   doubleGaussianB  ( in.getUntrackedParameter<int>    ("doubleGaussianB", 1));
  int   useBtag          ( in.getUntrackedParameter<int>    ("useBtag",         0));
  int   selectByBTagShape( in.getUntrackedParameter<int>    ("selectByBTagShape",0));
  int   doTypeBTag4      ( in.getUntrackedParameter<int>    ("doTypeBTag4", 0));
  int   doTypeBTag5      ( in.getUntrackedParameter<int>    ("doTypeBTag5", 0));
  int   doTypeBTag6      ( in.getUntrackedParameter<int>    ("doTypeBTag6", 0));
  int   doType0          ( in.getUntrackedParameter<int>    ("doType0", 0));
  int   doType1          ( in.getUntrackedParameter<int>    ("doType1", 0));
  int   doType2          ( in.getUntrackedParameter<int>    ("doType2", 0));
  int   doType3          ( in.getUntrackedParameter<int>    ("doType3", 0));
  int   doType6          ( in.getUntrackedParameter<int>    ("doType6", 0));
  int   doType7          ( in.getUntrackedParameter<int>    ("doType7", 0));
  int   doType0ByBTagShape(in.getUntrackedParameter<int>    ("doType0ByBTagShape", 0));
  int   doType1ByBTagShape(in.getUntrackedParameter<int>    ("doType1ByBTagShape", 0));
  int   doType2ByBTagShape(in.getUntrackedParameter<int>    ("doType2ByBTagShape", 0));
  int   doType3ByBTagShape(in.getUntrackedParameter<int>    ("doType3ByBTagShape", 0));
  int   doType6ByBTagShape(in.getUntrackedParameter<int>    ("doType6ByBTagShape", 0));
  int   useME            ( in.getParameter<int>             ("useME")     );
  int   useJac           ( in.getParameter<int>             ("useJac")    );
  int   useMET           ( in.getParameter<int>             ("useMET")    );
  int   useTF            ( in.getParameter<int>             ("useTF")     );
  int   usePDF           ( in.getParameter<int>             ("usePDF")    );
  int   norm             ( in.getUntrackedParameter<int>    ("norm",      0));
  int   hypo             ( in.getUntrackedParameter<int>    ("hypo",      0));
  int   SoB              ( in.getUntrackedParameter<int>    ("SoB",       1));
  int   doJERbias        ( in.getUntrackedParameter<int>    ("doJERbias", 0));

  int   doCSVup          ( in.getUntrackedParameter<int>    ("doCSVup",   0));
  int   doCSVdown        ( in.getUntrackedParameter<int>    ("doCSVdown", 0));
  int   doJECup          ( in.getUntrackedParameter<int>    ("doJECup",   0));
  int   doJECdown        ( in.getUntrackedParameter<int>    ("doJECdown", 0));
  int   doJERup          ( in.getUntrackedParameter<int>    ("doJERup",   0));
  int   doJERdown        ( in.getUntrackedParameter<int>    ("doJERdown", 0));
  int   fixNumEvJob      ( in.getUntrackedParameter<int>    ("fixNumEvJob",1));
  vector<string> functions(in.getParameter<vector<string> > ("functions"));
  vector<int>    evLimits (in.getParameter<vector<int> >    ("evLimits"));

  int   print            ( in.getParameter<int>             ("printout")   );
  int   debug            ( in.getParameter<int>             ("debug")      );


  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@ INITIALIZE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

  // flag to discriminate data/MC
  bool isMC = true;

  // upper and lower bounds to be processed
  int evLow  = evLimits[0];
  int evHigh = evLimits[1];

  TStopwatch* clock = new TStopwatch();

  // file with b-tag pdf
  TFile* fCP = TFile::Open(pathToCP.c_str(),"READ");

  // b-tag pdf for b-quark ('b'), c-quark ('c'), and light jets ('l')
  map<string,TH1F*> btagger;   if( useBtag && fCP!=0 ){
    btagger["b_Bin0"] = fCP->Get("csv_b_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_b_Bin0__csvReco") : 0;
    btagger["b_Bin1"] = fCP->Get("csv_b_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_b_Bin1__csvReco") : 0;
    btagger["c_Bin0"] = fCP->Get("csv_c_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_c_Bin0__csvReco") : 0;
    btagger["c_Bin1"] = fCP->Get("csv_c_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_c_Bin1__csvReco") : 0;
    btagger["l_Bin0"] = fCP->Get("csv_l_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_l_Bin0__csvReco") : 0;
    btagger["l_Bin1"] = fCP->Get("csv_l_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_l_Bin1__csvReco") : 0;
  }
  else if( useBtag && fCP!=0 ){
    cout << "Cound not find " << pathToCP << ": exit" << endl;
    return 0;
  }

  // Higgs mass values for scan
  const int nHiggsMassPoints  = massesH.size();
  double mH[nHiggsMassPoints];
  for( unsigned int m = 0; m < massesH.size() ; m++)
    mH[m] = massesH[m];

  // Top mass values for scan
  const int nTopMassPoints  = massesT.size();
  double mT[nTopMassPoints]; 
  for( unsigned int m = 0; m < massesT.size() ; m++)
    mT[m] = massesT[m];


  // not supported...
  if( nHiggsMassPoints>1 && nTopMassPoints>1){
    cout << "Cannot handle two mass scans at the same time... return." << endl;
    return 1;
  }
  if( nHiggsMassPoints>NMAXMASS || nTopMassPoints>NMAXMASS){
    cout << "Too many mass points required... return" << endl;
    return 1;
  }


  // configure MEIntegrator
  MEIntegratorNew* meIntegrator = new MEIntegratorNew( pathToTF , 4 , int(verbose));
  if( norm == 0)
    meIntegrator->setWeightNorm( MEIntegratorNew::None );
  else if( norm==1 )
    meIntegrator->setWeightNorm( MEIntegratorNew::xSec );
  else if( norm==2 )
    meIntegrator->setWeightNorm( MEIntegratorNew::Acc );
  else{
    cout << "Unsupported normalization... exit" << endl;
    delete meIntegrator;
    return 0;
  }

  // set normalization formulas ( not used if norm==0 )
  meIntegrator->setNormFormulas( TString(functions[0].c_str()),  
				 TString(functions[1].c_str()),  
				 TString(functions[2].c_str()),
				 TString(functions[3].c_str()),  
				 TString(functions[4].c_str()),  
				 TString(functions[5].c_str()),
				 TString(functions[6].c_str())
				 );
  
  // initialize top and W mass
  meIntegrator->setTopMass( MT , MW );
  meIntegrator->setTopSign( 1 );
  //  meIntegrator->setTopSign( leptCharge);

  // configure ME calculation
  meIntegrator->setUseME (useME);
  meIntegrator->setUseJac(useJac);
  meIntegrator->setUseMET(useMET);
  meIntegrator->setUseTF (useTF);
  meIntegrator->setUsePDF(usePDF);

  // use double-gaussian for b quark energy TF
  meIntegrator->setUseRefinedTF(doubleGaussianB);

  // use nominal TF from ROOT file
  meIntegrator->initTFparameters(1.0,1.0,1.0,1.0, 1.0);
  if(switchoffOL){
    meIntegrator->switchOffOL(); 
    cout << "*** Switching off OpenLoops to speed-up the calculation ***" << endl;
  }



  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@ OUTPUT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


  // clean output file (if any)
  gSystem->Exec(("rm "+outFileName).c_str());



  // output file
  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");

  // total event counter for normalization
  TH1F*  hcounter = new TH1F("hcounter","",1,0,1);
  int events_     = 0;
  
  // output tree
  TTree* tree  = new TTree("tree","");

  // counts how many events have been analyzed (to have fixed-size jobs)
  int counter_;
  // number of jet-quark permutations (for signal and bkg hypotheses)
  int nPermut_, nPermut_alt_;
  // number of (Higgs/Top) mass points
  int nMassPoints_;
  // total number integration per event (nPermut*nMassPoints)
  int nTotInteg_,nTotInteg_alt_;
  // number of matches to higgs quarks among tagged jets
  int matchesH_;
  // number of matches to W quarks among un-tagged jets
  int matchesW_;
  // number of matches to top quarks among tagged jets
  int matchesT_;
  // number of matches to higgs quarks among all jets
  int matchesHAll_;
  // number of matches to W quarks among all jets
  int matchesWAll_;
  // number of matches to top quarks among all jets
  int matchesTAll_;
  // count how many quarks from W decay overlap by dR<0.5
  int overlapLight_;
  // count how many b-quarks overlap by dR<0.5
  int overlapHeavy_;
  // integration type
  int type_;
  // num. of b-hadrons and c-quarks
  int nSimBs_; //, nC_, nCTop_;
  // num of b-hadrons inside jets
  int nMatchSimBs_;
  // type-dependent flags
  int flag_type0_;
  int flag_type1_;
  int flag_type2_;
  int flag_type3_;
  int flag_type4_;
  int flag_type6_;

  // event-wise **ME** probability (summed over permutations)
  // ...for ttH...
  float probAtSgn_;
  // ...for ttbb...
  float probAtSgn_alt_;

  // event-wise **ME X btag** probability (summed over permutations)
  // ...for ttH...
  float probAtSgn_ttbb_;
  // ...for tt(bb,jj,bj,cc)...
  float probAtSgn_alt_ttbb_;
  float probAtSgn_alt_ttjj_;
  float probAtSgn_alt_ttbj_;
  float probAtSgn_alt_ttcc_;

  // per-permutation and mass value **ME** probability
  float probAtSgn_permut_       [NMAXPERMUT*NMAXMASS];
  float probAtSgn_alt_permut_   [NMAXPERMUT*NMAXMASS];
  float probAtSgnErr_permut_    [NMAXPERMUT*NMAXMASS];
  float probAtSgnErr_alt_permut_[NMAXPERMUT*NMAXMASS];

  // per-permutation **btag** probability
  float probAtSgn_bb_permut_[NMAXPERMUT];
  float probAtSgn_bj_permut_[NMAXPERMUT];
  float probAtSgn_cc_permut_[NMAXPERMUT];
  float probAtSgn_jj_permut_[NMAXPERMUT];

  // masses to be scanned
  float mH_scan_[NMAXMASS];
  float mT_scan_[NMAXMASS];

  // event-dependent weight (for normalization)
  float weight_;
  // cpu time
  float time_;
  // event information
  EventInfo EVENT_;
  // num of PVs
  int nPVs_;
  // pu reweighting
  float PUweight_, PUweightP_, PUweightM_;
  // number of gen jets
  float lheNj_;

  // lepton kinematic (at most two leptons)
  int nLep_;
  float lepton_pt_    [2];
  float lepton_eta_   [2];
  float lepton_phi_   [2];
  float lepton_m_     [2];
  float lepton_charge_[2];
  float lepton_rIso_  [2];

  // met kinematic
  float MET_pt_;
  float MET_phi_;
  float MET_sumEt_;

  // jet kinematics (as passed via **jets** collection)
  int nJet_;
  float jet_pt_  [NMAXJETS];
  float jet_eta_ [NMAXJETS];
  float jet_phi_ [NMAXJETS];
  float jet_m_   [NMAXJETS];
  float jet_csv_ [NMAXJETS];

  // number of selected jets passing CSV L,M,T
  int numBTagL_, numBTagM_, numBTagT_;
  // nummber of selected jets
  int numJets_;
  // btag likelihood ratio
  float btag_LR_;

  // permutation -> jets association
  int   perm_to_jet_    [NMAXPERMUT];
  int   perm_to_jet_alt_[NMAXPERMUT];
  // permutation -> gen association
  int   perm_to_gen_     [NMAXPERMUT];
  int   perm_to_gen_alt_ [NMAXPERMUT];

  tree->Branch("counter",      &counter_,       "counter/I");
  tree->Branch("nPermut_s",    &nPermut_,       "nPermut_s/I");
  tree->Branch("nPermut_b",    &nPermut_alt_,   "nPermut_b/I");
  tree->Branch("nMassPoints",  &nMassPoints_,   "nMassPoints/I");  
  tree->Branch("nTotInteg_s",  &nTotInteg_,     "nTotInteg_s/I");  
  tree->Branch("nTotInteg_b",  &nTotInteg_alt_, "nTotInteg_b/I");  
  tree->Branch("matchesH",     &matchesH_,      "matchesH/I");
  tree->Branch("matchesW",     &matchesW_,      "matchesW/I");
  tree->Branch("matchesT",     &matchesT_,      "matchesT/I");
  tree->Branch("matchesHAll",  &matchesHAll_,   "matchesHAll/I");
  tree->Branch("matchesWAll",  &matchesWAll_,   "matchesWAll/I");
  tree->Branch("matchesTAll",  &matchesTAll_,   "matchesTAll/I");
  tree->Branch("overlapLight", &overlapLight_,  "overlapLight/I");
  tree->Branch("overlapHeavy", &overlapHeavy_,  "overlapHeavy/I");
  tree->Branch("type",         &type_,          "type/I");
  tree->Branch("nSimBs",       &nSimBs_,        "nSimBs/I");
  tree->Branch("nMatchSimBs",  &nMatchSimBs_,   "nMatchSimBs/I");
  //tree->Branch("nC",           &nC_,            "nC/I");
  //tree->Branch("nCTop",        &nCTop_,         "nCTop/I");
  tree->Branch("weight",       &weight_,        "weight/F");
  tree->Branch("time",         &time_,          "time/F");
  tree->Branch("flag_type0",   &flag_type0_,    "flag_type0/I");
  tree->Branch("flag_type1",   &flag_type1_,    "flag_type1/I");
  tree->Branch("flag_type2",   &flag_type2_,    "flag_type2/I");
  tree->Branch("flag_type3",   &flag_type3_,    "flag_type3/I");
  tree->Branch("flag_type4",   &flag_type4_,    "flag_type4/I");
  tree->Branch("flag_type6",   &flag_type6_,    "flag_type6/I");
  tree->Branch("EVENT",        &EVENT_,         "run/I:lumi/I:event/I:json/I");
  tree->Branch("nPVs",         &nPVs_,          "nPVs/I");
  tree->Branch("PUweight",     &PUweight_,      "PUweight/F");
  tree->Branch("PUweightP",    &PUweightP_,     "PUweightP/F");
  tree->Branch("PUweightM",    &PUweightM_,     "PUweightM/F");
  tree->Branch("lheNj",        &lheNj_,         "lheNj/F");


  // marginalized over permutations (all <=> Sum_{p=0}^{nTotPermut}( mH=MH, mT=MT ) )
  tree->Branch(Form("p_%d_all_s",     int(MH)),   &probAtSgn_,           Form("p_%d_all_s/F",              int(MH)) );
  tree->Branch(Form("p_%d_all_b",     int(MH)),   &probAtSgn_alt_,       Form("p_%d_all_b/F",              int(MH)) );
  tree->Branch(Form("p_%d_all_s_ttbb",int(MH)),   &probAtSgn_ttbb_,      Form("p_%d_all_s_ttbb/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttbb",int(MH)),   &probAtSgn_alt_ttbb_,  Form("p_%d_all_b_ttbb/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttjj",int(MH)),   &probAtSgn_alt_ttjj_,  Form("p_%d_all_b_ttjj/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttbj",int(MH)),   &probAtSgn_alt_ttbj_,  Form("p_%d_all_b_ttbj/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttcc",int(MH)),   &probAtSgn_alt_ttcc_,  Form("p_%d_all_b_ttcc/F",         int(MH)) );

  // differential in permutations and mass value. E.g.:
  //  p_vsMH_s[0],          ..., p_vsMH_s[  nPermut-1] => prob. per-permutation and for mH = masses[0]
  //  p_vsMH_s[nPermut], ...,    p_vsMH_s[2*nPermut-1] => prob. per-permutation and for mH = masses[1]
  //  ....
  tree->Branch("p_vsMH_s",     probAtSgn_permut_,       "p_vsMH_s[nTotInteg_s]/F" );
  tree->Branch("p_vsMT_b",     probAtSgn_alt_permut_,   "p_vsMT_b[nTotInteg_b]/F");
  tree->Branch("p_vsMH_Err_s", probAtSgnErr_permut_,    "p_vsMH_Err_s[nTotInteg_s]/F" );
  tree->Branch("p_msMT_Err_b", probAtSgnErr_alt_permut_,"p_vsMT_Err_b[nTotInteg_b]/F" );

  // differential in permutation
  tree->Branch("p_tt_bb",      probAtSgn_bb_permut_,    "p_tt_bb[nPermut_s]/F");
  tree->Branch("p_tt_bj",      probAtSgn_bj_permut_,    "p_tt_bj[nPermut_b]/F");
  tree->Branch("p_tt_cc",      probAtSgn_cc_permut_,    "p_tt_cc[nPermut_b]/F");
  tree->Branch("p_tt_jj",      probAtSgn_jj_permut_,    "p_tt_jj[nPermut_b]/F");
  
  // # of mass points scanned
  tree->Branch("mH_scan",       mH_scan_,          "mH_scan[nMassPoints]/F");
  tree->Branch("mT_scan",       mT_scan_,          "mT_scan[nMassPoints]/F");
 
  // lepton kinematics
  tree->Branch("nLep",                    &nLep_,        "nLep/I");
  tree->Branch("lepton_pt",               lepton_pt_,    "lepton_pt[nLep]/F");
  tree->Branch("lepton_eta",              lepton_eta_,   "lepton_eta[nLep]/F");
  tree->Branch("lepton_phi",              lepton_phi_,   "lepton_phi[nLep]/F");
  tree->Branch("lepton_m",                lepton_m_,     "lepton_m[nLep]/F");
  tree->Branch("lepton_charge",           lepton_charge_,"lepton_charge[nLep]/F");
  tree->Branch("lepton_rIso",             lepton_rIso_,  "lepton_rIso[nLep]/F");

  // MET kinematics
  tree->Branch("MET_pt",                  &MET_pt_,    "MET_pt/F");
  tree->Branch("MET_phi",                 &MET_phi_,   "MET_phi/F");
  tree->Branch("MET_sumEt",               &MET_sumEt_, "MET_sumEt/F");

  // jet kinematics
  tree->Branch("nJet",                    &nJet_,      "nJet/I");
  tree->Branch("jet_pt",                  jet_pt_,     "jet_pt[nJet]/F");
  tree->Branch("jet_eta",                 jet_eta_,    "jet_eta[nJet]/F");
  tree->Branch("jet_phi",                 jet_phi_,    "jet_phi[nJet]/F");
  tree->Branch("jet_m",                   jet_m_,      "jet_m[nJet]/F");
  tree->Branch("jet_csv",                 jet_csv_,    "jet_csv[nJet]/F");

  // Jet multiplicity
  tree->Branch("numBTagL",                &numBTagL_,    "numBTagL/I");
  tree->Branch("numBTagM",                &numBTagM_,    "numBTagM/I");
  tree->Branch("numBTagT",                &numBTagT_,    "numBTagT/I");
  tree->Branch("numJets",                 &numJets_,     "numJets/I");
  tree->Branch("btag_LR",                 &btag_LR_,     "btag_LR/F");

  // a map that associates to each permutation [p=0...nTotPermut] to the corresponding jet,
  // indexed according to the order in the jet_* collection
  //  E.g.: perm_to_jets[0] = 234567 <==> the first permutation (element '0' of p_vsMH_s/p_vsMT_b )
  //        associates element '2' of jets_* to the bLep, element '3' to W1Had, '4' to 'W2Had', '5' to bHad, and '6','7'
  //        to non-top radiation
  tree->Branch("perm_to_jet_s",             perm_to_jet_,    "perm_to_jet[nPermut_s]/I");
  tree->Branch("perm_to_gen_s" ,            perm_to_gen_,    "perm_to_gen[nPermut_s]/I");
  tree->Branch("perm_to_jet_b",             perm_to_jet_alt_,"perm_to_jet[nPermut_b]/I");
  tree->Branch("perm_to_gen_b" ,            perm_to_gen_alt_,"perm_to_gen[nPermut_b]/I");


  //tree->Branch("",                        _,  "[]/F");

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@ OPEN FILES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
 

  // read input files  
  bool openAllFiles  = false;
  Samples* mySamples = new Samples(openAllFiles, pathToFile, ordering, samples, lumi, verbose);
  vector<string> mySampleFiles;

  if(mySamples->IsOk()){

    cout << "Ok!" << endl;
    mySampleFiles = mySamples->Files();

    for( unsigned int i = 0 ; i < mySampleFiles.size(); i++){
      string sampleName       = mySampleFiles[i];

      if(verbose){
	cout << mySampleFiles[i] << " ==> " << mySamples->GetXSec(sampleName) 
	     << " pb,"
	     << " ==> weight = "            << mySamples->GetWeight(sampleName) << endl;
      }
    }
  }
  else{
    cout << "Problems... leaving" << endl;
    return 0;
  }

  // loop over input files
  for(unsigned int sample = 0 ; sample < mySampleFiles.size(); sample++){
    
    string currentName       = mySampleFiles[sample];

    if(currentName.find("Data")!=string::npos) isMC = false;

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;
    float scaleFactor        = mySamples->GetWeight(currentName);

    // variables to be used from the input files
    genParticleInfo genB, genBbar;
    genTopInfo      genTop, genTbar;
    metInfo         METtype1p2corr;
    EventInfo       EVENT;
    int             nvlep, nSimBs, /*nC, nCTop,*/ nhJets, naJets, nPVs;
    float           PUweight, PUweightP, PUweightM;
    float           lheNj;
    Int_t   vLepton_type      [999];
    Float_t vLepton_mass      [999];
    Float_t vLepton_pt        [999];
    Float_t vLepton_eta       [999];
    Float_t vLepton_phi       [999];
    Float_t vLepton_charge    [999];
    Float_t vLepton_pfCorrIso [999];
    Float_t hJet_pt           [999];
    Float_t hJet_eta          [999];
    Float_t hJet_phi          [999];
    Float_t hJet_e            [999];
    Float_t hJet_puJetIdL     [999];
    Float_t hJet_csv_nominal  [999];
    Float_t hJet_csv_upBC     [999];
    Float_t hJet_csv_downBC   [999];
    Float_t hJet_csv_upL      [999];
    Float_t hJet_csv_downL    [999];
    Float_t hJet_JECUnc       [999];
    Float_t hJet_genPt        [999];
    Float_t hJet_genEta       [999];
    Float_t hJet_genPhi       [999];
    Float_t aJet_pt           [999];
    Float_t aJet_eta          [999];
    Float_t aJet_phi          [999];
    Float_t aJet_e            [999];
    Float_t aJet_puJetIdL     [999];
    Float_t aJet_csv_nominal  [999];
    Float_t aJet_csv_upBC     [999];
    Float_t aJet_csv_downBC   [999];
    Float_t aJet_csv_upL      [999];
    Float_t aJet_csv_downL    [999];
    Float_t aJet_JECUnc       [999];
    Float_t aJet_genPt        [999];
    Float_t aJet_genEta       [999];
    Float_t aJet_genPhi       [999];
    float SimBsmass           [999];
    float SimBspt             [999];
    float SimBseta            [999];
    float SimBsphi            [999];
    //int nSvs;
    //float SvmassSv [999];
    //float Svpt     [999];
    //float Sveta    [999];
    //float Svphi    [999];
  
    currentTree->SetBranchAddress("EVENT",            &EVENT);
    currentTree->SetBranchAddress("PUweight",         &PUweight);
    currentTree->SetBranchAddress("PUweightP",        &PUweightP);
    currentTree->SetBranchAddress("PUweightM",        &PUweightM);
    currentTree->SetBranchAddress("lheNj",            &lheNj);    
    currentTree->SetBranchAddress("nhJets",           &nhJets);
    currentTree->SetBranchAddress("naJets",           &naJets);
    currentTree->SetBranchAddress("nSimBs",           &nSimBs);
    currentTree->SetBranchAddress("nvlep",            &nvlep);
    currentTree->SetBranchAddress("nPVs",             &nPVs);
    currentTree->SetBranchAddress("genB",             &genB);
    currentTree->SetBranchAddress("genBbar",          &genBbar);
    currentTree->SetBranchAddress("genTop",           &genTop);
    currentTree->SetBranchAddress("genTbar",          &genTbar);
    currentTree->SetBranchAddress("METtype1p2corr",   &METtype1p2corr);
    currentTree->SetBranchAddress("vLepton_charge",   vLepton_charge);
    currentTree->SetBranchAddress("vLepton_mass"  ,   vLepton_mass);
    currentTree->SetBranchAddress("vLepton_pt"    ,   vLepton_pt);
    currentTree->SetBranchAddress("vLepton_eta"   ,   vLepton_eta);
    currentTree->SetBranchAddress("vLepton_phi"   ,   vLepton_phi);
    currentTree->SetBranchAddress("vLepton_charge",   vLepton_charge);
    currentTree->SetBranchAddress("vLepton_pfCorrIso",vLepton_pfCorrIso);
    currentTree->SetBranchAddress("vLepton_type",     vLepton_type);
    currentTree->SetBranchAddress("hJet_pt",          hJet_pt);    
    currentTree->SetBranchAddress("hJet_eta",         hJet_eta);    
    currentTree->SetBranchAddress("hJet_phi",         hJet_phi);    
    currentTree->SetBranchAddress("hJet_e",           hJet_e);    
    currentTree->SetBranchAddress("hJet_puJetIdL",    hJet_puJetIdL);
    currentTree->SetBranchAddress("hJet_csv_nominal", hJet_csv_nominal);
    currentTree->SetBranchAddress("hJet_csv_upBC",    hJet_csv_upBC);
    currentTree->SetBranchAddress("hJet_csv_downBC",  hJet_csv_downBC);
    currentTree->SetBranchAddress("hJet_csv_upL",     hJet_csv_upL);
    currentTree->SetBranchAddress("hJet_csv_downL",   hJet_csv_downL);
    currentTree->SetBranchAddress("hJet_JECUnc",      hJet_JECUnc);
    currentTree->SetBranchAddress("hJet_genPt",       hJet_genPt);
    currentTree->SetBranchAddress("hJet_genEta",      hJet_genEta);
    currentTree->SetBranchAddress("hJet_genPhi",      hJet_genPhi);
    currentTree->SetBranchAddress("aJet_pt",          aJet_pt);    
    currentTree->SetBranchAddress("aJet_eta",         aJet_eta);    
    currentTree->SetBranchAddress("aJet_phi",         aJet_phi);    
    currentTree->SetBranchAddress("aJet_e",           aJet_e);    
    currentTree->SetBranchAddress("aJet_puJetIdL",    aJet_puJetIdL);
    currentTree->SetBranchAddress("aJet_csv_nominal", aJet_csv_nominal);
    currentTree->SetBranchAddress("aJet_csv_upBC",    aJet_csv_upBC);
    currentTree->SetBranchAddress("aJet_csv_downBC",  aJet_csv_downBC);
    currentTree->SetBranchAddress("aJet_csv_upL",     aJet_csv_upL);
    currentTree->SetBranchAddress("aJet_csv_downL",   aJet_csv_downL);
    currentTree->SetBranchAddress("aJet_JECUnc",      aJet_JECUnc);
    currentTree->SetBranchAddress("aJet_genPt",       aJet_genPt);
    currentTree->SetBranchAddress("aJet_genEta",      aJet_genEta);
    currentTree->SetBranchAddress("aJet_genPhi",      aJet_genPhi);
    currentTree->SetBranchAddress("SimBs_mass",   SimBsmass);
    currentTree->SetBranchAddress("SimBs_pt",     SimBspt);
    currentTree->SetBranchAddress("SimBs_eta",    SimBseta);
    currentTree->SetBranchAddress("SimBs_phi",    SimBsphi);
    //currentTree->SetBranchAddress("nSvs",        &nSvs);
    //currentTree->SetBranchAddress("Sv_massSv",   SvmassSv);
    //currentTree->SetBranchAddress("Sv_pt",       Svpt);
    //currentTree->SetBranchAddress("Sv_eta",      Sveta);
    //currentTree->SetBranchAddress("Sv_phi",      Svphi);



 
    /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
    /* @@@@@@@@@@@@@@@@@@@@@@ FWLITE EVENT LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */



    genParton genLFquark, genHiggs;

    int counter = 0;
    int nentries = 0;

    for(unsigned int iFile = 0; iFile<inFileNames.size(); ++iFile){

      
      TFile* fwfile = TFile::Open(inFileNames[iFile].c_str());
   


      // TFile* fwfile = TFile::Open("./root/el_edmntuples.root");
      //TFile fwfile("./root/singleTopEdmNtuple_TChannel.root");
      if(fwfile){
      	cout<<"file opened"<<endl;
      }
      fwlite::Event evt(fwfile);

      //    for(unsigned int  i= 0; i<2; i++){
      
      
      // loop over entries

      

    for(evt.toBegin(); !evt.atEnd(); ++evt){
      nentries++;
      //!!! questo e' da spostare BUILDING PARTONS

      vector<JetObservable> jetsEDM;
      vector<TLorentzVector> tightMus, tightEls, vetoEls, vetoMus;
      bool properEventSL;

      //Set the following flag to test just one interesting event
      // if(counter==0) continue;
      // if(counter==1) break;

     
      // reset variables
      probAtSgn_          =  0.;
      probAtSgn_alt_      =  0.;
      probAtSgn_ttbb_     =  0.;
      probAtSgn_alt_ttbb_ =  0.;
      probAtSgn_alt_ttbj_ =  0.;
      probAtSgn_alt_ttcc_ =  0.;
      probAtSgn_alt_ttjj_ =  0.;

      nMassPoints_        = TMath::Max(nHiggsMassPoints,nTopMassPoints);
      flag_type0_         = -99;
      flag_type1_         = -99;
      flag_type2_         = -99;
      flag_type3_         = -99;
      flag_type4_         = -99;
      flag_type6_         = -99;

      btag_LR_            = -99;

      for( int k = 0; k < 2; k++){
	lepton_pt_[k] = -99; lepton_eta_[k] = -99; lepton_phi_[k] = -99; lepton_m_[k] = -99; lepton_charge_[k] = -99; lepton_rIso_[k] = -99;
      }
      for( int k = 0; k < NMAXJETS; k++){
	jet_pt_[k] = -99; jet_eta_[k] = -99; jet_phi_[k] = -99; jet_m_[k] = -99;  jet_csv_[k] = -99;
      }
      for( int k = 0; k < NMAXPERMUT; k++){
	perm_to_jet_    [k] = -99;
	perm_to_gen_    [k] = -99;
	perm_to_jet_alt_[k] = -99;
	perm_to_gen_alt_[k] = -99;
      }


      // save the values into the tree (save mH[0] in case no scan is required)
      for( unsigned int m = 0; m < (unsigned int)nHiggsMassPoints ; m++){
	mH_scan_[m] = mH[m];
      }
      for( unsigned int t = 0; t < (unsigned int)nTopMassPoints ; t++){
	mT_scan_[t] = mT[t];
      }
      
      

      //*************************************************************
      //--------------------- RECO JETS -----------------------------
      //*************************************************************


      
      fwlite::Handle<std::vector<float>> jPt;
      jPt.getByLabel(evt, "nTupleTopJetsPF","topJetsPFPt");
      fwlite::Handle<std::vector<float>> jE;
      jE.getByLabel(evt, "nTupleTopJetsPF","topJetsPFE");
      if(doJECup){
	fwlite::Handle<std::vector<float>> jPt;
	jPt.getByLabel(evt, "nTupleTopJetsPF","topJetsPFPtJESUp");
	fwlite::Handle<std::vector<float>> jE;
	jE.getByLabel(evt, "nTupleTopJetsPF","topJetsPFEJESUp");
      }
      else if(doJECdown){
	fwlite::Handle<std::vector<float>> jPt;
	jPt.getByLabel(evt, "nTupleTopJetsPF","topJetsPFPtJESDown");
	fwlite::Handle<std::vector<float>> jE;
	jE.getByLabel(evt, "nTupleTopJetsPF","topJetsPFEJESDown");
      }
      else if(doJERup){
	fwlite::Handle<std::vector<float>> jPt;
	jPt.getByLabel(evt, "nTupleTopJetsPF","topJetsPFPtJERUp");
	fwlite::Handle<std::vector<float>> jE;
	jE.getByLabel(evt, "nTupleTopJetsPF","topJetsPFEJERUp");
      }
      else if(doJERdown){
	fwlite::Handle<std::vector<float>> jPt;
	jPt.getByLabel(evt, "nTupleTopJetsPF","topJetsPFPtJERDown");
	fwlite::Handle<std::vector<float>> jE;
	jE.getByLabel(evt, "nTupleTopJetsPF","topJetsPFEJERDown");
      }

      fwlite::Handle<std::vector<float>> jEta;
      jEta.getByLabel(evt, "nTupleTopJetsPF","topJetsPFEta");
      fwlite::Handle<std::vector<float>> jPhi;
      jPhi.getByLabel(evt, "nTupleTopJetsPF","topJetsPFPhi");
      fwlite::Handle<std::vector<float>> jCSV;
      jCSV.getByLabel(evt, "nTupleTopJetsPF","topJetsPFCombinedSecondaryVertexBJetTags");


//       cout<<"size :"<<jPt->size()<<endl;
//       cout<<"pt :"<<(*jPt)[0]<<endl;
//       cout<<"Pt :"<<evt.getRun().id()<<endl;

      for(unsigned int k = 0; k<jPt->size(); ++k ){
	TLorentzVector p4;
	p4.SetPtEtaPhiE((*jPt)[k], (*jEta)[k], (*jPhi)[k], (*jE)[k]);
	JetObservable jet;
	jet.p4 = p4;
	jet.csv =  (*jCSV)[k];
	jetsEDM.push_back(jet);
      } // end loop over jet collection per event


      //*************************************************************
      //-------------------- RECO LEPTONS ---------------------------
      //*************************************************************

      // TIGHT LEPTONS
      fwlite::Handle<std::vector<float>> mPt;
      mPt.getByLabel(evt, "nTupleMuons","tightMuonsPt");
      fwlite::Handle<std::vector<float>> mEta;
      mEta.getByLabel(evt, "nTupleMuons","tightMuonsEta");
      fwlite::Handle<std::vector<float>> mPhi;
      mPhi.getByLabel(evt, "nTupleMuons","tightMuonsPhi");
      fwlite::Handle<std::vector<float>> mE;
      mE.getByLabel(evt, "nTupleMuons","tightMuonsE");
      fwlite::Handle<std::vector<float>> mCharge;
      mCharge.getByLabel(evt, "nTupleMuons","tightMuonsCharge");



      for(unsigned int k = 0; k<mPt->size(); ++k ){
	TLorentzVector p4;
	p4.SetPtEtaPhiE((*mPt)[k], (*mEta)[k], (*mPhi)[k], (*mE)[k]);
	tightMus.push_back(p4);
      } // end loop over muons collection per event

      int ntightMu = tightMus.size(); 


      fwlite::Handle<std::vector<float>> ePt;
      ePt.getByLabel(evt, "nTupleElectrons","tightElectronsPt");
      fwlite::Handle<std::vector<float>> eEta;
      eEta.getByLabel(evt, "nTupleElectrons","tightElectronsEta");
      fwlite::Handle<std::vector<float>> ePhi;
      ePhi.getByLabel(evt, "nTupleElectrons","tightElectronsPhi");
      fwlite::Handle<std::vector<float>> eE;
      eE.getByLabel(evt, "nTupleElectrons","tightElectronsE");
      fwlite::Handle<std::vector<float>> eCharge;
      eCharge.getByLabel(evt, "nTupleElectrons","tightElectronsCharge");

      for(unsigned int k = 0; k<ePt->size(); ++k ){
	TLorentzVector p4;
	p4.SetPtEtaPhiE((*ePt)[k], (*eEta)[k], (*ePhi)[k], (*eE)[k]);
	tightEls.push_back(p4);
      } // end loop over electrons collection per event

      int ntightEl = tightEls.size(); 


      // VETO LEPTONS
      fwlite::Handle<std::vector<float>> vmPt;
      vmPt.getByLabel(evt, "nTupleVetoMuons","vetoMuonsPt");
      fwlite::Handle<std::vector<float>> vmEta;
      vmEta.getByLabel(evt, "nTupleVetoMuons","vetoMuonsEta");
      fwlite::Handle<std::vector<float>> vmPhi;
      vmPhi.getByLabel(evt, "nTupleVetoMuons","vetoMuonsPhi");
      fwlite::Handle<std::vector<float>> vmE;
      vmE.getByLabel(evt, "nTupleVetoMuons","vetoMuonsE");

      for(unsigned int k = 0; k<vmPt->size(); ++k ){
	TLorentzVector p4;
	p4.SetPtEtaPhiE((*vmPt)[k], (*vmEta)[k], (*vmPhi)[k], (*vmE)[k]);
	vetoMus.push_back(p4);
      } // end loop over veto muons collection per event

      int nvetoMu = vetoMus.size(); 

      fwlite::Handle<std::vector<float>> vePt;
      vePt.getByLabel(evt, "nTupleVetoElectrons","vetoElectronsPt");
      fwlite::Handle<std::vector<float>> veEta;
      veEta.getByLabel(evt, "nTupleVetoElectrons","vetoElectronsEta");
      fwlite::Handle<std::vector<float>> vePhi;
      vePhi.getByLabel(evt, "nTupleVetoElectrons","vetoElectronsPhi");
      fwlite::Handle<std::vector<float>> veE;
      veE.getByLabel(evt, "nTupleVetoElectrons","vetoElectronsE");

      for(unsigned int k = 0; k<vePt->size(); ++k ){
	TLorentzVector p4;
	p4.SetPtEtaPhiE((*vePt)[k], (*veEta)[k], (*vePhi)[k], (*veE)[k]);
	vetoEls.push_back(p4);
      } // end loop over veto electrons collection per event


      int nvetoEl = vetoEls.size(); 


      //*************************************************************
      //-------------------- RECO MET ---------------------------
      //*************************************************************

      fwlite::Handle<std::vector<float>> metPt;
      metPt.getByLabel(evt, "nTuplePatMETsPF","patMETsPFPt");
      fwlite::Handle<std::vector<float>> metPhi;
      metPhi.getByLabel(evt, "nTuplePatMETsPF","patMETsPFPhi");


      //*************************************************************
      //-------------------- GEN B FROM HIGGS -----------------------
      //*************************************************************

      fwlite::Handle<std::vector<float>> genBHPt;
      genBHPt.getByLabel(evt, "singleTopMCHiggsBQuark","MCHiggsBQuarkPt");
      fwlite::Handle<std::vector<float>> genBHEta;
      genBHEta.getByLabel(evt, "singleTopMCHiggsBQuark","MCHiggsBQuarkEta");
      fwlite::Handle<std::vector<float>> genBHPhi;
      genBHPhi.getByLabel(evt, "singleTopMCHiggsBQuark","MCHiggsBQuarkPhi");
      fwlite::Handle<std::vector<float>> genBHE;
      genBHE.getByLabel(evt, "singleTopMCHiggsBQuark","MCHiggsBQuarkE");
      fwlite::Handle<std::vector<float>> genBHCharge;
      genBHCharge.getByLabel(evt, "singleTopMCHiggsBQuark","MCHiggsBQuarkCharge");


      for(unsigned int k = 0; k<genBHPt->size(); ++k ){
	TLorentzVector p4;
	p4.SetPtEtaPhiE((*genBHPt)[k], (*genBHEta)[k], (*genBHPhi)[k], (*genBHE)[k]);
      } // end loop over gen B from Higgs collection per event


      //*************************************************************
      //-------------------- GEN B FROM TOP -------------------------
      //*************************************************************

      fwlite::Handle<std::vector<float>> genBTopPt;
      genBTopPt.getByLabel(evt, "singleTopMCTopsBQuark","MCtopsBQuarkPt");
      fwlite::Handle<std::vector<float>> genBTopEta;
      genBTopEta.getByLabel(evt, "singleTopMCTopsBQuark","MCtopsBQuarkEta");
      fwlite::Handle<std::vector<float>> genBTopPhi;
      genBTopPhi.getByLabel(evt, "singleTopMCTopsBQuark","MCtopsBQuarkPhi");
      fwlite::Handle<std::vector<float>> genBTopE;
      genBTopE.getByLabel(evt, "singleTopMCTopsBQuark","MCtopsBQuarkE");
      fwlite::Handle<std::vector<float>> genBTopCharge;
      genBTopCharge.getByLabel(evt, "singleTopMCTopsBQuark","MCtopsBQuarkCharge");


      for(unsigned int k = 0; k<genBTopPt->size(); ++k ){
	TLorentzVector p4;
	p4.SetPtEtaPhiE((*genBTopPt)[k], (*genBTopEta)[k], (*genBTopPhi)[k], (*genBTopE)[k]);
      } // end loop over gen B from Higgs collection per event



      //*************************************************************
      //-------------------- GEN LF QUARK ---------------------------
      //*************************************************************
      // ATT: REQUEST ON THE MOTHER!!!!


      fwlite::Handle<std::vector<float>> genLFqPt;
      genLFqPt.getByLabel(evt, "singleTopMCQuarks","MCquarksPt");
      fwlite::Handle<std::vector<float>> genLFqEta;
      genLFqEta.getByLabel(evt, "singleTopMCQuarks","MCquarksEta");
      fwlite::Handle<std::vector<float>> genLFqPhi;
      genLFqPhi.getByLabel(evt, "singleTopMCQuarks","MCquarksPhi");
      fwlite::Handle<std::vector<float>> genLFqE;
      genLFqE.getByLabel(evt, "singleTopMCQuarks","MCquarksE");





      for(unsigned int k = 0; k<genLFqPt->size(); ++k ){
	TLorentzVector p4;
	p4.SetPtEtaPhiE((*genLFqPt)[k], (*genLFqEta)[k], (*genLFqPhi)[k], (*genLFqE)[k]);


	genLFquark.mass = p4.M();
	genLFquark.pt = p4.Pt();
	genLFquark.phi = p4.Phi();
	genLFquark.eta = p4.Eta();
	genLFquark.momid1 = 0.;
	genLFquark.momid2 = 0.;
	genLFquark.status = 0.;
	genLFquark.charge = 0.;
	
      } // end loop over gen B from Higgs collection per event


      //*************************************************************
      //-------------------- GEN LEPTONS ---------------------------
      //*************************************************************

      fwlite::Handle<std::vector<float>> genLTopPt;
      genLTopPt.getByLabel(evt, "singleTopMCTopsLepton","MCtopsLeptonPt");
      fwlite::Handle<std::vector<float>> genLTopEta;
      genLTopEta.getByLabel(evt, "singleTopMCTopsLepton","MCtopsLeptonEta");
      fwlite::Handle<std::vector<float>> genLTopPhi;
      genLTopPhi.getByLabel(evt, "singleTopMCTopsLepton","MCtopsLeptonPhi");
      fwlite::Handle<std::vector<float>> genLTopE;
      genLTopE.getByLabel(evt, "singleTopMCTopsLepton","MCtopsLeptonE");
      fwlite::Handle<std::vector<float>> genLTopCharge;
      genLTopCharge.getByLabel(evt, "singleTopMCTopsLepton","MCtopsLeptonCharge");


      for(unsigned int k = 0; k<genLTopPt->size(); ++k ){
	TLorentzVector p4;
	p4.SetPtEtaPhiE((*genLTopPt)[k], (*genLTopEta)[k], (*genLTopPhi)[k], (*genLTopE)[k]);

      } // end loop over gen lepton from Top collection per event


      //*************************************************************
      //-------------------- GEN Neutrino ---------------------------
      //*************************************************************

      fwlite::Handle<std::vector<float>> genNuTopPt;
      genNuTopPt.getByLabel(evt, "singleTopMCTopsNeutrino","MCtopsNeutrinoPt");
      fwlite::Handle<std::vector<float>> genNuTopEta;
      genNuTopEta.getByLabel(evt, "singleTopMCTopsNeutrino","MCtopsNeutrinoEta");
      fwlite::Handle<std::vector<float>> genNuTopPhi;
      genNuTopPhi.getByLabel(evt, "singleTopMCTopsNeutrino","MCtopsNeutrinoPhi");
      fwlite::Handle<std::vector<float>> genNuTopE;
      genNuTopE.getByLabel(evt, "singleTopMCTopsNeutrino","MCtopsNeutrinoE");




      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@ LEPTON SELECTION @@@@@@@@@@@@@@@@@@@@@@@@@  */
      

      // charged leptons and MET
      TLorentzVector leptonLV;
      TLorentzVector neutrinoLV;


      ///////////////////////////////////
      //         SL events             //
      ///////////////////////////////////

      properEventSL = false;      
      if( ntightMu==1 && ntightEl == 0 && nvetoMu == 1 && nvetoEl == 0 ){

	// first lepton...
	leptonLV = tightMus[0] ;
	// cut on muons already applied at skim level
	properEventSL = true;

	// save lepton kinematics...
	nLep_ = 1;

	lepton_pt_     [0] = leptonLV.Pt();
	lepton_eta_    [0] = leptonLV.Eta();
	lepton_phi_    [0] = leptonLV.Phi();
	lepton_m_      [0] = leptonLV.M();
	lepton_charge_ [0] = (*mCharge)[0];

      }
      else if( ntightEl==1 && ntightMu == 0 && nvetoMu == 0 && nvetoEl == 1 ){
	
	leptonLV = tightEls[0] ;
	properEventSL = true;

	nLep_ = 1;

	lepton_pt_     [0] = leptonLV.Pt();
	lepton_eta_    [0] = leptonLV.Eta();
	lepton_phi_    [0] = leptonLV.Phi();
	lepton_m_      [0] = leptonLV.M();
	lepton_charge_ [0] = (*eCharge)[0];

      }


      ////////////////////////////////////////////////////////////////////////

      // MET
      float metPt_ = (*metPt)[0];
      float metPhi_ = (*metPhi)[1];
      float nuPx = metPt_*TMath::Cos(metPhi_);
      float nuPy = metPt_*TMath::Sin(metPhi_);
      float nuE  = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);
      neutrinoLV.SetPxPyPzE(nuPx,nuPy,0. ,nuE);

      // save MET kinematics into the tree...
      MET_pt_    = neutrinoLV.Pt();
      MET_phi_   = neutrinoLV.Phi();
      MET_sumEt_ = TMath::Sqrt(metPt_*metPt_);

      // continue if leptons do not satisfy cuts
      if( !(properEventSL) ) continue;



      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@ JET SELECTION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

      // this container will hold the jets
      std::vector<JetObservable> jet_map;

      int nJets = jetsEDM.size();
      // loop over jets
      for(int j = 0; j < nJets; j++){

	  float ptGen = -99.;
	  //	  if(coll==0 && hJet_genPt[hj]>0.) ptGen = hJet_genPt[hj];

	  float pt     = jetsEDM[j].p4.Pt();
	  float eta    = jetsEDM[j].p4.Eta();
	  float phi    = jetsEDM[j].p4.Phi();
	  float e      = jetsEDM[j].p4.E();
	  float m2     = e*e - pt*pt*TMath::CosH(eta)*TMath::CosH(eta);
	  if(m2<0) m2 = 0.; 
	  float m      = TMath::Sqrt( m2 ); 
	  float csv =  jetsEDM[j].csv;

	  // only jets in acceptance...
	  if( TMath::Abs(eta)> 4.7 ) continue;

	  // only jets above pt cut...
	  if( pt < 30  ) continue;	  

	  TLorentzVector p4;
	  p4.SetPtEtaPhiM( pt, eta, phi, m );

	  //!!! Missing btag SF uncertainties in edm ntuples.

// 	  // for csv systematics
// 	  float csv_nominal =  (coll==0) ? hJet_csv_nominal[hj] : aJet_csv_nominal[hj];
// 	  float csv_upBC    =  (coll==0) ? hJet_csv_upBC   [hj] : aJet_csv_upBC   [hj];
// 	  float csv_downBC  =  (coll==0) ? hJet_csv_downBC [hj] : aJet_csv_downBC [hj];
// 	  float csv_upL     =  (coll==0) ? hJet_csv_upL    [hj] : aJet_csv_upL    [hj];
// 	  float csv_downL   =  (coll==0) ? hJet_csv_downL  [hj] : aJet_csv_downL  [hj];
// 	  float csv = csv_nominal;
// 	  if     ( doCSVup  ) csv =  TMath::Max(csv_upBC,   csv_upL);
// 	  else if( doCSVdown) csv =  TMath::Min(csv_downBC, csv_downL);
// 	  else{}

	  // the b-tagger output 
	  //  ==> Min needed because in csvUp, csv can exceed 1..., 
	  //      Max needed because we crunch  [-inf,0[ -> {0.}

	  csv         =  TMath::Min( TMath::Max( csv, float(0.)), float(0.999999) );

	  // the jet observables (p4 and csv)
	  JetObservable myJet;
	  myJet.p4  = p4;
	  myJet.csv = csv; 
	  
	  // use pt to order jet collection
	  //jet_map[ p4.Pt() ] = myJet;
	  jet_map.push_back( myJet );

	  if( debug )
	  cout << "Jet #" << j << " => (" << pt << "," << eta << "," << phi << "," << m << "), ID=" << endl;
	      
      }// end loop over jets
      
      /* @@@@@@@@@@@@@@@@@@@@@@@@ JET ORDERING @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      
      // order jet list by Pt
      std::sort( jet_map.begin(), jet_map.end(), JetObservableListerByPt() );
      
      // if use btag shape, order by decreasing CSV 
      //       <=> when considering only a subset of the jets, this ensures that the combination
      //           obtained from CSVM only jets is among those considered
      if( selectByBTagShape ) 
	std::sort( jet_map.begin(), jet_map.end(), JetObservableListerByCSV() );

   
      // fill arrays of jets
      std::vector<TLorentzVector>  jets_p4;
      std::vector<double>          jets_csv;
      //std::vector<double>          jets_csv_prob_b;
      //std::vector<double>          jets_csv_prob_c;
      //std::vector<double>          jets_csv_prob_j;
      int jetsAboveCut = 0;

      //for( jet_map_it = jet_map.begin() ; jet_map_it != jet_map.end(); jet_map_it++){
      for(unsigned int jj = 0; jj < jet_map.size() ; jj++ ){
	
	// the four-vector
	TLorentzVector p4 = jet_map[jj].p4; //(jet_map_it->second).p4;
	
	// the csv value
	float csv = jet_map[jj].csv;
	
	// count jets above 40 GeV
	if( p4.Pt()>30 ) jetsAboveCut++;
	
	// store jet p4...
	jets_p4.push_back ( p4  );
	// store csv
	jets_csv.push_back( csv );
	
	// needed to find appropriate csv PDF
	string bin = "";
	if( TMath::Abs( p4.Eta() ) <= 1.0 ) 
	  bin = "Bin0";
	if( TMath::Abs( p4.Eta() ) >  1.0 ) 
	  bin = "Bin1";

	// store PDF(csv)
	//jets_csv_prob_b.push_back( btagger["b_"+bin]!=0 ? btagger["b_"+bin]->GetBinContent( btagger["b_"+bin]->FindBin( csv ) ) : 1.);
	//jets_csv_prob_c.push_back( btagger["c_"+bin]!=0 ? btagger["c_"+bin]->GetBinContent( btagger["c_"+bin]->FindBin( csv ) ) : 1.);
	//jets_csv_prob_j.push_back( btagger["l_"+bin]!=0 ? btagger["l_"+bin]->GetBinContent( btagger["l_"+bin]->FindBin( csv ) ) : 1.);

      }

      // continue if not enough jets
      if( jetsAboveCut<4 ){
	//cout << "Less then 4 jets.. continue" << endl;
	continue;
      }
      

      /* @@@@@@@@@@@@@@@@@@@@@@@@ JET ORDERING @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

      // jet multiplicity
      int numJets30UntagM = 0; 
      int numJets30UntagL = 0; 
      int numJets30BtagL  = 0; 
      int numJets30BtagM  = 0; 
      int numJets30BtagT  = 0;

      // this vector contains the indices of up to 6 jets passing the CSVM(L) b-tag selection... 
      // (type 7 uses two working points)
      vector<unsigned int> btag_indices;
      vector<unsigned int> btagLoose_indices;

      // this vector contains the indices of up to 6 jets failing the CSVM(L) b-tag selection... 
      vector<unsigned int> buntag_indices;
      vector<unsigned int> buntagLoose_indices;

      // this vector contains the indices of up to 12 jets passing ANY b-tag selection... 
      vector<unsigned int> banytag_indices;  
    

      for(unsigned int k = 0; k < jets_p4.size(); k++){  

	float csv_k = jets_csv[k];
	float pt_k  = (jets_p4[k]).Pt();

	// passes CSVL...
	int btag_L = csv_k>0.244 ;	
	// passes CSVM...
	int btag_M = csv_k>0.679 ;
	// passes CSVT...
	int btag_T = csv_k>0.898 ;	

	// any btag value...
	if( pt_k>30 ) banytag_indices.push_back( k );

	// passing at least CSVL...
	if( pt_k>30 &&  btag_L ){
	  numJets30BtagL ++;
	  btagLoose_indices.push_back( k );
	}
	// passing at least CSVM...
	if( pt_k>30 &&  btag_M ){
	  numJets30BtagM ++;
	  btag_indices.push_back( k );
	}
	// passing at least CSVT...
	if( pt_k>30 &&  btag_T ){
	  numJets30BtagT ++;
	}

	// failing at most CSVL...
	if( pt_k>30 && !btag_L ){
	  numJets30UntagL++;
	  buntagLoose_indices.push_back( k );
	}
	// failing at most CSVM...
	if( pt_k>30 && !btag_M ){
	  numJets30UntagM++;
	  buntag_indices.push_back( k );
	}

      }

      numBTagL_ = numJets30BtagL;
      numBTagM_ = numJets30BtagM;
      numBTagT_ = numJets30BtagT;
      numJets_  = ( numJets30BtagM + numJets30UntagM);


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ EVENT SELECTION @@@@@@@@@@@@@@@@@@@@@@@@@  */


      ////////////////////////////////////////////////////
      // SEMILEPTONIC EVENTS                            //
      ////////////////////////////////////////////////////

      // categories defined by jet and btagged jet multiplicity (Nj,Nb)
      bool analyze_th_4j3b     = properEventSL && numJets30BtagM==3 && numJets30UntagM==1  && doType0;
      bool analyze_th_AtLeast_4j3b     = properEventSL && numJets30BtagM==3 && numJets30UntagM>=1  && doType0;



      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ANALYSIS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


      //  input 4-vectors
      vector<TLorentzVector> jets;

      // internal map: [ position in "jets" ] -> [ position in "jets_p4" ]
      std::map< unsigned int, unsigned int> pos_to_index;

      // consider th event only if of the desired type
//       if( analyze_typeBTag6  || analyze_typeBTag5  || analyze_typeBTag4  || 
// 	  analyze_type0      || analyze_type1      || analyze_type2      || analyze_type3      || analyze_type6       || analyze_type7 ||
// 	  analyze_type0_BTag || analyze_type1_BTag || analyze_type2_BTag || analyze_type3_BTag || analyze_type6_BTag){	

      unsigned int indq = 5;

   //    cout<<"Event to be analyzed: "<<analyze_th_4j3b <<endl;

//       cout<<"Event to be analyzed: "<<analyze_th_AtLeast_4j3b<<endl;
 
      if( analyze_th_4j3b ){
	//      if( analyze_th_4j3b || analyze_th_AtLeast_4j3b){
	
	events_++;

	indq = 5;
	// find out which two untagged jets come from W->qq'

	//cout<<"MEAnalysisNew - The event has to be analyzed"<<endl;
	//cout<<"MEAnalysisNew - Print is set to: "<<print<<endl;
	//cout<<"MEAnalysisNew - Print fixnumber of ev per job: "<<fixNumEvJob<<endl;
	//cout<<"MEAnalysisNew - Print event counter : "<<(counter>=evLow && counter<=evHigh)<<endl;
	
	
	counter++;      
	//cout<<"MEAnalysisNew - Output of the condition to continue: "<<(fixNumEvJob && !(counter>=evLow && counter<=evHigh) )<<endl;
	
	
	//	if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	if(print) cout << "MEAnalysisNew - Processing event # " << counter << " (TYPE 4j 3b)." << " Event ID " << EVENT.event << endl;	   
	
	/////////////////////////////////////////////////////
	type_       =  0;
	nPermut_    =  3;
	nPermut_alt_=  6; //!!!to change	  
	meIntegrator->setIntType( MEIntegratorNew::SL4j3b );
	/////////////////////////////////////////////////////	  

      }
      
      else{ 
	//cout << "MEAnalysisNew - Inconsistency in the analysis... continue." << endl; 
	continue; 
      }    


	// DEBUG
      if(debug){
	  cout << "*** Event ID " << EVENT.event << " *** " << endl;
	  cout << " ==> SL=" << int(properEventSL) << endl;//", DL=" << properEventDL << endl;
	  cout << "     NJets " << numJets_ << " (" << numBTagM_ << " tagged)" << endl;
	  cout << "     b-tagged: " << endl;
	//   for( unsigned int jj = 0; jj<btag_indices_backup.size(); jj++)
// 	    cout << "     (" 
// 		 << jets_p4[ btag_indices_backup[jj] ].Pt() << "," 
// 		 << jets_p4[ btag_indices_backup[jj] ].Eta() << ","
// 		 << jets_p4[ btag_indices_backup[jj] ].Phi() << "," 
// 		 << jets_p4[ btag_indices_backup[jj] ].M() << "), CSV= " 
// 		 << jets_csv[btag_indices_backup[jj] ] << endl;
	  cout << "     b-untagged: " << endl;
// 	  for( unsigned int jj = 0; jj<buntag_indices_backup.size(); jj++)
// 	    cout << "     (" 
// 		 << jets_p4[ buntag_indices_backup[jj] ].Pt() << "," 
// 		 << jets_p4[ buntag_indices_backup[jj] ].Eta() << ","
// 		 << jets_p4[ buntag_indices_backup[jj] ].Phi() << "," 
// 		 << jets_p4[ buntag_indices_backup[jj] ].M() << "), CSV= " 
// 		 << jets_csv[buntag_indices_backup[jj] ] << endl;
	  cout << "     btag probability is " << btag_LR_ << endl;
	//   if(passes_btagshape){
// 	    cout << "     @@@@@ the jet collection has been re-ordered according to btag probability @@@@@@" << endl;
// 	    cout << "     b-tagged: " << endl;
// 	    for( unsigned int jj = 0; jj < btag_indices.size(); jj++)
// 	      cout << "     (" << jets_p4[ btag_indices[jj] ].Pt() << "," << jets_p4[ btag_indices[jj] ].Eta() << ","
// 		   << jets_p4[ btag_indices[jj] ].Phi() << "," << jets_p4[ btag_indices[jj] ].M() << "), CSV= " << jets_csv[ btag_indices[jj] ] << endl;
// 	    cout << "     b-untagged: " << endl;
// 	    for( unsigned int jj = 0; jj<buntag_indices.size(); jj++)
// 	      cout << "     (" << jets_p4[ buntag_indices[jj] ].Pt() << "," << jets_p4[ buntag_indices[jj] ].Eta() << ","
// 		   << jets_p4[ buntag_indices[jj] ].Phi() << "," << jets_p4[ buntag_indices[jj] ].M() << "), CSV= " << jets_csv[buntag_indices[jj] ] << endl;
// 	  }
	  
	  
	}
	
	

	// total number of integrations
	nTotInteg_      = nPermut_    * nMassPoints_;
	nTotInteg_alt_  = nPermut_alt_* nMassPoints_;
	
	
	// setup jet collection
	jets.clear();
	jets.push_back( leptonLV     );  
	jets.push_back( neutrinoLV   );  

	// b1,...,w1,w2 are indices for jets_p4 collection;
	// This is a map between the internal ordering. For  th  bTop=2, b1Higgs=3, b2Higgs = 4 and lfq = 5
	pos_to_index.clear();
	 

	//	cout<< "MEAnalysisNew - TYPE "<<type_<<endl; 
	if( type_==0){	  
	  jets.push_back( jets_p4[ btag_indices[0] ]);
	  jets.push_back( jets_p4[ btag_indices[1] ]);
	  jets.push_back( jets_p4[ btag_indices[2] ]);
	  jets.push_back( jets_p4[ buntag_indices[0] ]);
	  
	  pos_to_index[2] = btag_indices[0];
	  pos_to_index[3] = btag_indices[1]; // dummy
	  pos_to_index[4] = btag_indices[2]; // dummy
	  pos_to_index[5] = buntag_indices[0];
	}
	else{ /* ... */ }
	



	// save jet kinematics into the tree...
	nJet_ = jets.size();

	//cout<<"MEAnalysisNew - nPermut_ "<< nPermut_<<endl;
	//cout<<"MEAnalysisNew - nMassPoints_ "<< nMassPoints_<<endl;
	//cout<<"MEAnalysisNew - nTotInteg_ "<< nTotInteg_<<endl;
	//cout<<"MEAnalysisNew - nTotInteg_alt_ "<< nTotInteg_alt_<<endl;


	//cout<<"MEAnalysisNew - njets_ "<< nJet_<<endl;

	//cout<<"MEAnalysisNew - Filling jets"<<endl;


	for(int q = 0; q < nJet_ ; q++ ){
	  // kinematics
	  jet_pt_ [q] = jets[q].Pt(); 
	  jet_eta_[q] = jets[q].Eta(); 
	  jet_phi_[q] = jets[q].Phi(); 	    
	  jet_m_  [q] = jets[q].M(); 
	  jet_csv_[q] = q>1 ? jets_csv[ pos_to_index[q] ] : -99.;
	}
	

	//cout<<"MEAnalysisNew - Filled jets"<<endl;


	// set all prob. to 0.0;
	for(int p = 0 ; p < nTotInteg_; p++){
	  probAtSgn_permut_       [p] = 0.;
	  probAtSgnErr_permut_    [p] = 0.;
	}


	for(int p = 0 ; p < nPermut_; p++){
	  probAtSgn_bb_permut_ [p] = 0.;
	}


	//	nTotInteg_alt_ = 2;

	for(int p = 0 ; p < nTotInteg_alt_; p++){
	  probAtSgn_alt_permut_   [p] = 0.;
	  probAtSgnErr_alt_permut_[p] = 0.;
	}

	for(int p = 0 ; p < nPermut_alt_; p++){
	  probAtSgn_bj_permut_ [p] = 0.;
	  probAtSgn_cc_permut_ [p] = 0.;
	  probAtSgn_jj_permut_ [p] = 0.;
	}
	
	cout<<"initialization permutations done"<<endl;

	/////////////////////////////////////////////////////////////
	
	// check if there is a tag-untag pair that satisfies the "cs-tag" 
	for( unsigned int w = 0; w<btag_indices.size(); w++){
	  
	  float m1 = ( jets_p4[btag_indices[w]] + jets_p4[indq] ).M();
	  
	  if( (m1>(MwL+5) && m1<(MwH-5))  && type_== 0 ){
	    flag_type0_ = 0; 
	    if( jets_csv[btag_indices[w]]<0.95 ) flag_type0_ = 1; 
	    if( jets_csv[btag_indices[w]]<0.90 ) flag_type0_ = 2; 
	    if( jets_csv[btag_indices[w]]<0.85 ) flag_type0_ = 3;
	    if( jets_csv[btag_indices[w]]<0.80 ) flag_type0_ = 4;  
	  }
	}
	
	
	/////////////////////////////////////////////////////////////
	
	/////////////////////////////////////////////////////////////
	
	//cout<<"MEAnalysisNew - Set jets "<<endl;


	// init reco particles
	meIntegrator->setJets(&jets);	

	// init MET stuff
	meIntegrator->setSumEt(MET_sumEt_);
	meIntegrator->setMEtCov(-99,-99,0);
	// set top charge
	meIntegrator->setTopSign(lepton_charge_[0]);
	
	// specify if topLep has pdgid +6 or -6
	//meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
	  
	// start the clock...	  
	clock->Start();
	

	//cout<<"MEAnalysisNew - Starting loop over higgs masses"<<endl;

	// loop over Higgs mass values...
	for(int m = 0; m < nHiggsMassPoints ; m++){
	  meIntegrator->setMass( mH[m] );
	  
	  // loop over Top mass values...
	  for(int t = 0; t < nTopMassPoints ; t++){
	    meIntegrator->setTopMass( mT[t] , MW );
	    
	    // these are used for bookkeeping
	    double maxP_s = 0.;
	    double maxP_b = 0.;
	    
	    // number of permutations for which p has been calculated
	    int num_s = 0;
	    int num_b = 0;

	    // loop over hypothesis [ TH, TTJ ]
	    
	    for(int hyp = 0 ; hyp<2;  hyp++){


	      //cout<<"MEAnalysisNew - Starting loop over hypo"<<endl;	    


	      // choose which permutations to consider;
	      int* permutList = 0;
	      if     ( type_ == 0 ) permutList = hyp==0 ? permutations_th4j3b     : permutations_ttq4j3b;
	      //	      if( type_ == -2 ) permutList = hyp==0 ?  permutations_6J_S      : permutations_6J_B;
	      else{ 
		cout << "No permutations found...continue." << endl; 
		continue; 
	      }
	      

	      // loop over permutations
	      for(unsigned int pos = 0; pos < (unsigned int)( hyp==0 ? nPermut_ : nPermut_alt_ ) ; pos++){
		
		// consider permutation #pos & save permutation-to-jet mas into the tree...
		meIntegrator->initVersors( permutList[pos] );
		if( hyp==0 ) perm_to_jet_    [pos] =  permutList[pos];
		if( hyp==1 ) perm_to_jet_alt_[pos] =  permutList[pos];
		
		// index of the four jets associated to b-quarks or W->qq
		// for the bkg hypothesis, b1 and b2 are respectively the b quark from the hadronic top 
		// and the first quark from the  hadronic W of the haronic top

 
		int bLep_pos = (permutList[pos])%10000/1000;
		int b1_pos   = (permutList[pos])%1000/100;
		int b2_pos   = (permutList[pos])%100/10;       
		int lfq_pos  = (permutList[pos])%10/1;       
	
		//!!! MC Matching : to be restored
	
		// find if this particular permutation matches the expectation
	// 	int bLep_match = 0;
// 		if(TMath::Abs(TOPLEPB.Py())>0 && deltaR( jets_p4[ pos_to_index[bLep_pos] ], TOPLEPB ) < GENJETDR){
// 		  bLep_match = 1;
// 		}
// 		int w1_match = 0;  
// 		int w2_match = 0;
// 		if( (TMath::Abs(TOPHADW1.Py())>0 && deltaR( jets_p4[ pos_to_index[w1_pos] ], TOPHADW1 ) < GENJETDR) || 
// 		    (TMath::Abs(TOPHADW2.Py())>0 && deltaR( jets_p4[ pos_to_index[w1_pos] ], TOPHADW2 ) < GENJETDR)){
// 		  w1_match = 1;
// 		}
// 		if( (TMath::Abs(TOPHADW1.Py())>0 && deltaR( jets_p4[ pos_to_index[w2_pos] ], TOPHADW1 ) < GENJETDR) || 
// 		    (TMath::Abs(TOPHADW2.Py())>0 && deltaR( jets_p4[ pos_to_index[w2_pos] ], TOPHADW2 ) < GENJETDR)){
// 		  w2_match = 1;
// 		}
// 		int bHad_match = 0;
// 		if(TMath::Abs(TOPHADB.Py())>0 && deltaR( jets_p4[ pos_to_index[bHad_pos] ], TOPHADB ) < GENJETDR){
// 		  bHad_match = 1;
// 		}
// 		int b1_match = 0;  
// 		int b2_match = 0;
// 		if( (TMath::Abs(HIGGSB1.Py())>0 && deltaR( jets_p4[ pos_to_index[b1_pos] ], HIGGSB1 ) < GENJETDR) || 
// 		    (TMath::Abs(HIGGSB2.Py())>0 && deltaR( jets_p4[ pos_to_index[b1_pos] ], HIGGSB2 ) < GENJETDR)){
// 		  b1_match = 1;
// 		}
// 		if( (TMath::Abs(HIGGSB1.Py())>0 && deltaR( jets_p4[ pos_to_index[b2_pos] ], HIGGSB1 ) < GENJETDR) || 
// 		    (TMath::Abs(HIGGSB2.Py())>0 && deltaR( jets_p4[ pos_to_index[b2_pos] ], HIGGSB2 ) < GENJETDR)){
// 		  b2_match = 1;
// 		}
		
		// save an integer with this convention:
		//  > 1st digit = 1/0 if bLep candidate is correct/wrong
		//  > 2nd digit = 1/0 if matched/not matched to W-quark
		//  > ... 

		//		if( hyp==0 ) perm_to_gen_    [pos] = 100000*bLep_match + 10000*w1_match + 1000*w2_match + 100*bHad_match + 10*b1_match + 1*b2_match;
		//		if( hyp==1 ) perm_to_gen_alt_[pos] = 100000*bLep_match + 10000*w1_match + 1000*w2_match + 100*bHad_match + 10*b1_match + 1*b2_match;
		
		
		// check invariant mass of jet system:
		double mass, massLow, massHigh;
		bool skip        = !( meIntegrator->compatibilityCheck    (0.95, /*print*/ 0, mass, massLow, massHigh ) );
		//bool skip_WHad   = false;
		//bool skip_TopHad = false;
		//		if( type_==0 || type_==3 ){
		//skip_WHad   = !( meIntegrator->compatibilityCheck_WHad  (0.98, /*print*/ 0, mass, massLow, massHigh ) );
		//skip_TopHad = !( meIntegrator->compatibilityCheck_TopHad(0.98, /*print*/ 0, mass, massLow, massHigh ) );
		//}
		
		
		// if use btag, determine the b-tag probability density
	// 	if( useBtag ) {	       			
// 		  double p_b_bLep =  jets_csv_prob_b[ pos_to_index[bLep_pos] ];
// 		  double p_b_bHad =  jets_csv_prob_b[ pos_to_index[bHad_pos] ];
// 		  double p_b_b1   =  jets_csv_prob_b[ pos_to_index[b1_pos] ];
// 		  double p_c_b1   =  jets_csv_prob_c[ pos_to_index[b1_pos] ];
// 		  double p_j_b1   =  jets_csv_prob_j[ pos_to_index[b1_pos] ];
// 		  double p_b_b2   =  jets_csv_prob_b[ pos_to_index[b2_pos] ];
// 		  double p_c_b2   =  jets_csv_prob_c[ pos_to_index[b2_pos] ];
// 		  double p_j_b2   =  jets_csv_prob_j[ pos_to_index[b2_pos] ];
// 		  double p_j_w1   =  1.0;
// 		  double p_j_w2   =  1.0;
		  
// 		  // the b-tag probability for the two untagged jets...
// 		  if( type_==0 || type_==3 || type_==-2){
// 		    p_j_w1 = jets_csv_prob_j[ pos_to_index[w1_pos] ];
// 		    p_j_w2 = jets_csv_prob_j[ pos_to_index[w2_pos] ];
// 		  }
// 		  // the b-tag probability for the one untagged jet...
// 		  else if( type_==1 || type_==2 || type_==-1){
// 		    p_j_w1 = jets_csv_prob_j[ pos_to_index[w1_pos] ];
// 		    p_j_w2 = 1.0;
// 		  }
// 		  // there are untagged jets...
// 		  else{
// 		    p_j_w1 = 1.0;
// 		    p_j_w2 = 1.0;
// 		  }
		  
// 		  // DEBUG
// 		  if( debug && 0 ){
// 		    cout << "Hyp=" << hyp <<"  [BTag M="   << numBTagM_ << "]" << endl;
// 		    cout << " bLep: p_b("   << jets_csv[pos_to_index[bLep_pos]] << ")=" << p_b_bLep;
// 		    cout << " bHad: p_b("   << jets_csv[pos_to_index[bHad_pos]] << ")=" << p_b_bHad;
// 		    cout << " b1  : p_b("   << jets_csv[pos_to_index[b1_pos]]   << ")=" << p_b_b1;
// 		    cout << " b2  : p_b("   << jets_csv[pos_to_index[b2_pos]]   << ")=" << p_b_b2;
// 		    //		    if(!(type_>=6 || type_==-3)) 
// 		    //cout << " w1  : p_j(" << jets_csv[pos_to_index[w1_pos]]   << ")=" << p_j_w1;
// 		    //if(!(type_>=6 || type_==-3 || type_==1 || type_==2)) 
// 		    //cout << " w2  : p_j(" << jets_csv[pos_to_index[w2_pos]]   << ")=" << p_j_w2 << endl;
// 		    cout << " P = "         <<  p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 * p_j_w1 * p_j_w2 << endl;
// 		  }
		  
// 		  // fill arrays with per-permutation probability
// 		  if(hyp==0){
// 		    probAtSgn_bb_permut_[pos] =  p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 * p_j_w1 * p_j_w2;
// 		  }
// 		  if(hyp==1){
// 		    probAtSgn_bj_permut_[pos] =  p_b_bLep * p_b_bHad * (p_b_b1 * p_j_b2 + p_j_b1 * p_b_b2 )*0.5 * p_j_w1 * p_j_w2;
// 		    probAtSgn_cc_permut_[pos] =  p_b_bLep * p_b_bHad * p_c_b1 * p_c_b2 * p_j_w1 * p_j_w2;
// 		    probAtSgn_jj_permut_[pos] =  p_b_bLep * p_b_bHad * p_j_b1 * p_j_b2 * p_j_w1 * p_j_w2;		  
// 		  }

// 		}
		
		// if doing scan over b-tag only, don't need amplitude...
		if( type_<0  ) continue;
		
		
		// if type 0/3 and incompatible with MW or MT (and we are not scanning vs MT) continue
		// ( this applies to both hypotheses )
		// if( nTopMassPoints==1 && (skip_WHad || skip_TopHad) ){    
// 		  if(print) cout << "Skip (THad check failed)...          Perm. #" << pos << endl;
// 		  continue;
// 		}
		
	
		
		//cout<<"MEAnalysisNew - Setting Ranges"<<endl;
	
		// retrieve intergration boundaries from meIntegrator
		//pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
		//cout<<"MEAnalysisNew - First range defind "<<endl;
		pair<double, double> range_x1 =  make_pair(-1,1);
		//cout<<"MEAnalysisNew - Second range defind "<<endl;
		//pair<double, double> range_x2 =  make_pair(-PI,PI);
		//cout<<"MEAnalysisNew - Third range defind "<<endl;	    
		pair<double, double> range_x3 =  make_pair(-1,1);
		//cout<<"MEAnalysisNew - Fourth range defind "<<endl;
		pair<double, double> range_x4 =  useMET ? (meIntegrator->getNuPhiCI(0.95)) : make_pair(-PI,PI);
		//cout<<"MEAnalysisNew - Fifth range defind "<<endl;
		pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
		//cout<<"MEAnalysisNew - Sisth range defind "<<endl;
		//pair<double, double> range_x6 = (meIntegrator->getB2EnergyCI(0.95));
		//cout<<"MEAnalysisNew - Eighth range defind "<<endl;
		pair<double, double> range_x7 = (meIntegrator->getLFQEnergyCI(0.95));
		//cout<<"MEAnalysisNew - Nineth range defind "<<endl;
		pair<double, double> range_x8 = make_pair(-PI,PI);

		// boundaries
		//double x0L = range_x0.first; double x0U = range_x0.second;
		double x1L = range_x1.first; double x1U = range_x1.second;
		//double x2L = range_x2.first; double x2U = range_x2.second;
		double x3L = range_x3.first; double x3U = range_x3.second;
		double x4L = range_x4.first; double x4U = range_x4.second;
		double x5L = range_x5.first; double x5U = range_x5.second;
		//double x6L = range_x6.first; double x6U = range_x6.second;
		double x7L = range_x7.first; double x7U = range_x7.second;
		double x8L = range_x8.first; double x8U = range_x8.second;
	      
		// these hold for the sgn integration and type0...
		// 		double xLmode0_s[4] = {x0L, x3L, x4L, x5L};
		// 		double xUmode0_s[4] = {x0U, x3U, x4U, x5U};


		// new type 0: 4j3b Integrate over Eq, Eb(higgs), Phi and CosTheta nu
		double xLmode0_s[4] = {x7L, x3L, x4L, x5L};
		double xUmode0_s[4] = {x7U, x3U, x4U, x5U};
		// these hold for the bkg integration and type0...
		double xLmode0_b[6] = { x7L, x1L, x4L, x5L, x3L, x8L };
		double xUmode0_b[6] = { x7U, x1U, x4U, x5U, x3U, x8U};		 


		//double xLmode0_b[5] = {x0L, x3L, x4L, x5L, x6L};
		//double xUmode0_b[5] = {x0U, x3U, x4U, x5U, x6U};		 
	      

		//cout<<"MEAnalysisNew - Ranges Set "<<endl;


		// number of integration variables (TTH hypothesis)
		int nParam;
		if     ( type_==0 )  nParam = 4; // Eq, eta_nu, phi_nu, Eb
		else if( type_==1 )  nParam = 6; // Eq, eta_q', phi_q', eta_nu, phi_nu, Eb
		else if( type_==2 )  nParam = 6; // Eq, eta_nu, phi_nu, Eb,eta_q', phi_q'
		else if( type_==3 )  nParam = 4; // Eq, eta_nu, phi_nu, Eb
		else{
		  cout << "No type match...continue." << endl;
		  continue;
		}
		
		// the per-permutation probability...
		double p    = 0.;	
		// the per-permutation probability error...
		double pErr = 0.;	       	
		
	     	  	
		// if doing higgs mass scan, don't consider bkg hypo
		if( nHiggsMassPoints>1 && hyp==1 ) continue;
		
		// if doing top mass scan, don't consider sgn hypo
		if( nTopMassPoints>1 && hyp==0 ) continue;
		
		// if consider only one hypothesis (SoB=0) 
		// and the current hypo is not the desired one, continue...
		if( SoB==0 && hyp!=hypo) continue;
		
		// if current hypo is TH, but M(b1b2) incompatible with 125
		// (and we are not scanning vs MH) continue...
		if( hyp==0 && nHiggsMassPoints==1 && skip){
		  if(print){
		    cout << "Skip    hypo " << (hyp==0 ? "tH " : "ttj") 
			 << " [MH=" << mH[m] << ", MT=" << mT[t]
			 << "] Perm. #" << pos;
		    cout << " => p=" << p << endl;
		  }
		  continue;
		}
		
		// increment counters
		if(hyp==0) num_s++;
		if(hyp==1) num_b++;

		// if NOT doing top mass scan, and hyp==1
		// and no permutations for hyp==0 accepted, continue (the weight will be 0 anyway)
	

		if( nTopMassPoints==1 && hyp==1 && num_s==0){
		  if(print){
		    cout << "Skip    hypo " << (hyp==0 ? "tH " : "ttj") 
			 << " because no valid permutations found" << endl;
		  }
		  continue;
		}

		// setup hypothesis
		if(print){
		  cout << "Testing hypo " << (hyp==0 ? "tH " : "ttj") 
		       << " [MH=" << mH[m] << ", MT=" << mT[t]
		       << "] Perm. #" << pos;	       
		}
		//cout<<"MEAnalysisNew - setting hypo"<<endl;
		meIntegrator->setHypo(hyp);
		//cout<<"MEAnalysisNew - hypo set"<<endl;
		// initial number of function calles
		int intPoints = 4000;
		if( type_==0 )  intPoints =  2000;
		//		if( type_==0 )  intPoints =  2000;
		if( type_==1 )  intPoints =  4000;
		if( type_==2 )  intPoints =  4000;
		if( type_==3 )  intPoints =  2000;
		//		if( type_==6 )  intPoints = 10000;
		//		if( type_==7 )  intPoints = 10000;
		
		  // count how many time the integration is rerun per permutation
		int ntries = 0;
		
		// skip ME calculation... for debugging
		if(speedup==0){
		  
		  // refinement: redo integration if it returned a bad chi2
		  while( ntries < MAX_REEVAL_TRIES){
		    
		    // integrand
		    ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, nParam+(2*hyp));
		    // VEGAS integrator
		    ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		    ig2.SetFunction(toIntegrate);
		    // setup # of parameters
		    meIntegrator->SetPar( nParam+(2*hyp) );	 
		    
		    // the integration ranges depend on hyp and type
		    if     ( type_==0 ){
		      //cout<<"MEAnalysisNew - Integral"<<endl;
		      p = (hyp==0 ? ig2.Integral(xLmode0_s, xUmode0_s) : ig2.Integral(xLmode0_b, xUmode0_b));
		    }
		    else{ /* ... */ }		    
		    
		    // chi2 of the integration
		    double chi2 =  ig2.ChiSqr();
		    
		    // error from the last VEGAS iteration
		    pErr =  ig2.Error();
		    
		    // check if the actual permutation returned a small or large number...
		    // if the value is less than 10% of the largest found that far, skip 
		    // the improvement
		    if( hyp==0 ){
		      if( p>maxP_s ) maxP_s = p;		  
		      else if( p<0.1*maxP_s) ntries = MAX_REEVAL_TRIES+1;
		    }
		    else{
		      if( p>maxP_b ) maxP_b = p;		  
		      else if( p<0.1*maxP_b) ntries = MAX_REEVAL_TRIES+1;	
		    }
		    
		    // if the chi2 is bad, increse # of points and repeat the integration
		    if( chi2 > maxChi2_ ){
		      ntries++;
		      intPoints *= 1.5;
		    }
		    // otherwise, just go to the next permutation...
		    else 
		      ntries = MAX_REEVAL_TRIES+1;			
		  }	
		}
		else{
		  // can still be interested in b-tagging, so set p=1...
		  p = 1.;
		}
		  
		// to avoid problems
		if( TMath::IsNaN(p) )    p    = 0.;
		if( TMath::IsNaN(pErr) ) pErr = 0.;
		if(print) cout << " => p=(" << p << " +/- " << pErr << ")" << endl;
		
		/////////////////////////////////////////////////////////////
		  
		if( useBtag ){
		  
		  // total ME*btag probability for nominal MH and MT, sgn hypo
		  if( hyp==0 && mH[m]<MH+1.5 && mH[m]>MH-1.5 && mT[t]<MT+1.5 && mT[t]>MT-1.5){
		    probAtSgn_ttbb_         += ( p * probAtSgn_bb_permut_[pos] );
		  }
		  
		  // total ME*btag probability for nominal MH and MT, bkg hypo
		  if( hyp==1 && mH[m]<MH+1.5 && mH[m]>MH-1.5 && mT[t]<MT+1.5 && mT[t]>MT-1.5){
		    probAtSgn_alt_ttbb_     += ( p * probAtSgn_bb_permut_[pos] );
		    probAtSgn_alt_ttbj_     += ( p * probAtSgn_bj_permut_[pos] );
		    probAtSgn_alt_ttcc_     += ( p * probAtSgn_cc_permut_[pos] );
		    probAtSgn_alt_ttjj_     += ( p * probAtSgn_jj_permut_[pos] );		      
		  }		      		      
		}
		
		/////////////////////////////////////////////////////////////
		
		// save TH prob. per permutation AND mH
		if( hyp==0 && mT[t]<MT+1.5 && mT[t]>MT-1.5){
		  probAtSgn_permut_   [(unsigned int)(pos+m*nPermut_)] = p;
		  probAtSgnErr_permut_[(unsigned int)(pos+m*nPermut_)] = pErr;
		}
		// save TTJ prob. per permutation AND mT
		if( hyp==1 && mH[m]<MH+1.5 && mH[m]>MH-1.5 ){
		  probAtSgn_alt_permut_   [(unsigned int)(pos+t*nPermut_alt_)] = p;
		  probAtSgnErr_alt_permut_[(unsigned int)(pos+t*nPermut_alt_)] = pErr;
		}
		
		// total and per-permutation ME probability for nominal MH and MT
		if( mH[m]<MH+1.5 && mH[m]>MH-1.5 && mT[t]<MT+1.5 && mT[t]>MT-1.5){
		  if(hyp==0){
		    probAtSgn_     += p;
		  }
		  else{
		    probAtSgn_alt_ += p;
		  }
		}
		
		/////////////////////////////////////////////////////////////
		  
	      }  // permutations		  
	    }  // hypotheses
	  }  // nTopMassPoints
	}  // nHiggsMassPoints	    
	
	// stop clock and reset
	clock->Stop();
	time_ = clock->RealTime();
	clock->Reset();	    
	
	// this time is integrated over all hypotheses, permutations, and mass scan values
	if(print) cout << "Done in " << time_ << " sec" << endl;
	//}
    
      ///////////////////////////////////////////////////
      // ALL THE REST...                               //
      ///////////////////////////////////////////////////

	else{
	  continue;
	}
   
        counter_     = counter;
	nSimBs_      = nSimBs;
	weight_      = scaleFactor;
	nPVs_        = nPVs;
	
	EVENT_.run   = EVENT.run;
	EVENT_.lumi  = EVENT.lumi;
	EVENT_.event = EVENT.event;
	EVENT_.json  = EVENT.json;
	
// 	PUweight_    = PUweight;
// 	PUweightP_   = PUweightP;
// 	PUweightM_   = PUweightM;
	
// 	lheNj_       = lheNj;
      
	// fill the tree...
	tree->Fill();          



	
    } // end fwlite loop (on edmntuple)



  } // end loop over iFile

    // this histogram keeps track of the fraction of analyzed events per sample
    hcounter->SetBinContent(1,float(events_)/nentries);    
    
  } // samples


  // save the tree and the counting histo in the ROOT file
  fout_tmp->cd();

  hcounter->Write("",TObject::kOverwrite );
  tree->Write("",TObject::kOverwrite );
  fout_tmp->Close();

  // delete MEIntegrator and other allocated objects
  cout << "Delete meIntegrator..." << endl;
  delete meIntegrator;
  delete clock;
  // delete probAtSgn_permut_;
  //delete probAtSgnErr_permut_;
  //delete probAtSgn_alt_permut_;
  //delete probAtSgnErr_permut_;
  cout << "Finished!!!" << endl;
  
  // Done!!!
  return 0;
}
