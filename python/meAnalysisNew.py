import FWCore.ParameterSet.Types  as CfgTypes
import FWCore.ParameterSet.Config as cms

## ########################################

import FWCore.ParameterSet.VarParsing as opts


options = opts.VarParsing ('analysis')
del options._register['outputFile']
del options._singletons['outputFile']
import sys
print sys.argv

options.register ('outputLabel',
                  '',
                  opts.VarParsing.multiplicity.singleton,
                  opts.VarParsing.varType.string,
                  'Output label')

options.register ('sample',
                  '',
                  opts.VarParsing.multiplicity.singleton,
                  opts.VarParsing.varType.string,
                  'Sample analyzed')

options.register ('maxFiles',
                  0, # default value
                  opts.VarParsing.multiplicity.singleton, # singleton or list
                  opts.VarParsing.varType.int,          # string, int, or float
                  'Number of files process')

options.register ('skipFiles',
                  0, # default value
                  opts.VarParsing.multiplicity.singleton, # singleton or list
                  opts.VarParsing.varType.int,          # string, int, or float
                  'Number of files to skip')

options.register ('verbose',
                  False, # default value
                  opts.VarParsing.multiplicity.singleton, # singleton or list
                  opts.VarParsing.varType.bool,          # string, int, or float
                  'Enable/Disable verbosity')


options.parseArguments()

#if options.outputLabel: options.outputLabel += '_'

if options.skipFiles != 0:
    options._lists['inputFiles'] = options.inputFiles[options.skipFiles:]
   
if options.maxFiles != 0:
    n = min(options.maxFiles,len(options.inputFiles))
    options._lists['inputFiles'] = options.inputFiles[:n]
    

print '-'*80
print '- InputFiles'
print '-'*80
print '\n'.join(options.inputFiles)
print '-'*80

## ########################################



VType     = "_VType2"

xsecTT_SL = 103.0
xsecTT_FL = 24.8

process = cms.Process("MEAnalysisNew")

process.fwliteInput = cms.PSet(

    inFileNames   = cms.vstring(options.inputFiles),
    outFileName   = cms.string("/scratch/decosa/tH/MEM/"+options.sample+"/mem_"+options.outputLabel+".root"),
    pathToTF      = cms.string("./root/transferFunctionsNew_partonE_new.root"),
    pathToCP      = cms.string("./root/ControlPlotsNew_new.root"),
    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt"+VType+"/v2/"),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),

    samples       = cms.VPSet(


    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTH_HToBB_M-120_8TeV-pythia6'+VType),
    nickName = cms.string('TTH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1470*0.648)
    ),

    
    
    ),


    #SLNoBLep 10000
    #SLNoBHad 10000
    #SLNoHiggs 8000
    #SL2wj     2000
    #SL1wj     4000
    #DL       10000
    
    hypo          = cms.untracked.int32(0),
    SoB           = cms.untracked.int32(1),
    maxChi2       = cms.untracked.double(2.5), # 2.5
    norm          = cms.untracked.int32(0),

     functions     = cms.vstring(
    '8.95351e+18*TMath::Landau(x, 5.67600e+01,1.01258e+01)',                # incl
    '2.98474e+17*TMath::Landau(x,7.40196e+01 ,1.80142e+01)',                # SL4j3b ### to redefine, just keeping SL1wj
    '2.95547e+17*TMath::Landau(x,7.61581e+01 ,1.89245e+01)',                # SL2wj
    '2.98474e+17*TMath::Landau(x,7.40196e+01 ,1.80142e+01)',                # SL1wj
    '6.28300e+16*TMath::Landau(x,8.03060e+01 ,1.81679e+01)',                # SLNoBHad
    'x>150?2.44515e+27*x^(-5.35628e+00):1.24208e+18*exp((-3.63162e-02)*x)', # SLNoHiggs
    'x>=12 ? x^(-2.010e-01)*exp((-1.5785e-02)*x) : 4.184e-02*x'),           # tth Pt
 
    switchoffOL   = cms.untracked.int32(0), ###### CHECK HERE
    speedup       = cms.untracked.int32(0), ###### CHECK HERE

    doTypeBTag6   = cms.untracked.int32(0),  #SL 6 jets
    doTypeBTag5   = cms.untracked.int32(0),  #SL 5 jets
    doTypeBTag4   = cms.untracked.int32(0),  #DL 4 jets
    
    
    doType0       = cms.untracked.int32(1),  #SL(4,3) for tH
    doType1       = cms.untracked.int32(0),  #SL(4,2)  w/o W-tag
    doType2       = cms.untracked.int32(0),  #SL(4,1)
    doType3       = cms.untracked.int32(0),  #SL(4,3) 
    doType4       = cms.untracked.int32(0),  #SL(3,2)
    doType6       = cms.untracked.int32(0),  #DL(4,X)
    doType7       = cms.untracked.int32(0),  #DL(3M+1L,X)

    doType0ByBTagShape = cms.untracked.int32(0),
    doType1ByBTagShape = cms.untracked.int32(0),
    doType2ByBTagShape = cms.untracked.int32(0),
    doType3ByBTagShape = cms.untracked.int32(0),
    doType6ByBTagShape = cms.untracked.int32(0),

#    useME         = cms.int32(0),
    useME         = cms.int32(1),
    useJac        = cms.int32(1),
    useMET        = cms.int32(1),
    useTF         = cms.int32(1),
    usePDF        = cms.int32(1),


    doubleGaussianB  = cms.untracked.int32(1),
    useBtag          = cms.untracked.int32(0),
    selectByBTagShape= cms.untracked.int32(0),
    
    printout     = cms.int32(1),
    #    debug        = cms.int32(1),
    debug        = cms.int32(0),   
    verbose      = cms.bool(options.verbose),
#    verbose      = cms.bool(False),

    MH           = cms.untracked.double(125.00),
    MT           = cms.untracked.double(174.30),
    MW           = cms.untracked.double( 80.19),

    MwL          = cms.untracked.double(60),
    MwH          = cms.untracked.double(100),
    MhL          = cms.untracked.double(110),
    MhH          = cms.untracked.double(140),

    btag_prob_cut_6jets = cms.untracked.double(0.988),
    btag_prob_cut_5jets = cms.untracked.double(0.992),
    btag_prob_cut_4jets = cms.untracked.double(0.992),
    
    massesH      = cms.vdouble(125.),
   # massesH      = cms.vdouble(65., 75., 85., 95., 105., 115., 125., 135., 145., 155., 165., 185., 205., 225., 250., 300.),
    massesT      = cms.vdouble(174.3),
    #massesT      = cms.vdouble(145, 155, 165, 174, 185, 195, 205),

    fixNumEvJob    = cms.untracked.int32(0),
    evLimits       = cms.vint32(0,1),

    doJERbias  = cms.untracked.int32(0),   
    doCSVup    = cms.untracked.int32(0),
    doCSVdown  = cms.untracked.int32(0),
    doJECup    = cms.untracked.int32(0),
    doJECdown  = cms.untracked.int32(0),
    doJERup    = cms.untracked.int32(0),
    doJERdown  = cms.untracked.int32(0),

    )
