#!/bin/bash

FIRST=0
#LAST=10
LAST=141
STEP=1
#SAMPLE=tH_minus
SAMPLE=TTJets_SemiLeptMGDecays_8TeV-madgraph
LABEL=16Dec_


for SKIP in $( seq $FIRST $STEP $LAST)
do


  echo $SKIP

#  qexe_sh.py  -q all.q -t th_minus_$LABEL$SKIP -- "mkdir -p /scratch/decosa/tH/MEM/$SAMPLE; MEAnalysisNew ../python/meAnalysisNew.py inputFiles_load=th_minus_4FS.txt maxFiles=$STEP skipFiles=$SKIP outputLabel=$SKIP sample=$SAMPLE;  lcg-cp -b -D srmv2 /scratch/decosa/tH/MEM/$SAMPLE/mem_$SKIP.root  srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat/store/user/decosa/tH/MEM/$LABEL/$SAMPLE/mem_$SKIP.root" 

  qexe_sh.py   -q all.q -t ttjets_$LABEL$SKIP -- "mkdir -p /scratch/decosa/tH/MEM/$SAMPLE; MEAnalysisNew ../python/meAnalysisNew.py inputFiles_load=ttjets.txt maxFiles=$STEP skipFiles=$SKIP outputLabel=$SKIP sample=$SAMPLE;  lcg-cp -b -D srmv2 /scratch/decosa/tH/MEM/$SAMPLE/mem_$SKIP.root  srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat/store/user/decosa/tH/MEM/$LABEL/$SAMPLE/mem_$SKIP.root" 
  


#  qexe_sh.py -t th_minus_$SKIP -- MEAnalysisNew ../python/meAnalysisNew.py inputFiles_load=th_minus_4FS.txt maxFiles=$STEP skipFiles=$SKIP outputLabel=$SKIP 

done
