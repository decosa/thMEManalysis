import time
import numpy as np
start = time.time()

import sys
sys.argv.append('-b')
import ROOT
import commands, os

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleFont(42)
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)


c = ROOT.TCanvas()
ROOT.SetOwnership(c,False)
c.cd()
                

# Indicate the path on the t3 where to find mem root files
prefix = "dcap://t3se01.psi.ch:22125/"
path_s = "/pnfs/psi.ch/cms/trivcat/store/user/decosa/tH/MEM/13Dec_/tH_minus/"
path_b = "/pnfs/psi.ch/cms/trivcat/store/user/decosa/tH/MEM/16Dec_/TTJets_SemiLeptMGDecays_8TeV-madgraph/"
cmd_ = "ls -lrt " + path_s
cmd_s = "ls -lrt " + path_s
cmd_b = "ls -lrt " + path_b

# Set the legend
leg = ROOT.TLegend(.65, .7, .9, .85)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.SetTextFont(42)

# Pick up list of files
status_,ls_la_ = commands.getstatusoutput( cmd_ )
list_ = ls_la_.split(os.linesep)
list_.remove(list_[0])

files = []

for a in list_:
    b = a.split(" ")
    files.append(b[-1])

def getChain(path):
    # Pick up list of files
    cmd_ = "ls -lrt " + path
    status,ls_la_ = commands.getstatusoutput( cmd_ )
    list_ = ls_la_.split(os.linesep)
    list_.remove(list_[0])
    
    files = []
    
    for a in list_:
        b = a.split(" ")
        files.append(b[-1])

    chain = ROOT.TChain("tree")
    
    for f in files:
        chain.Add(prefix+path+f)
    return chain


def drawPs(chain):
    nEntries = chain.GetEntries()

    for i in xrange(nEntries):
        chain.GetEntry(i)
        ps = chain.p_125_all_s
        if(ps!=0 ):
            hs_.Fill(ps)
    return hs_

chain_s = getChain(path_s)
chain_b = getChain(path_b)

#hs_ = drawPs


hs_ = ROOT.TH1F("hypo S", "hypos", 100, 0.,1.)
hb_ = ROOT.TH1F("hypo B", "hypob", 100, 0.,1.)
hs_.SetFillColor(ROOT.kRed-3)
hs_.SetLineColor(ROOT.kRed+1)
hs_.SetFillStyle(3005)
#hb_.SetFillColor(ROOT.kGray)
#hb_.SetLineColor(ROOT.kGray+2)
hb_.SetFillColor(ROOT.kBlue-3)
hb_.SetLineColor(ROOT.kBlue+1)


ROOT.SetOwnership(hs_,False)
ROOT.SetOwnership(hb_,False)


nEntries = chain_s.GetEntries()
print nEntries

for i in xrange(nEntries):
    chain_s.GetEntry(i)
    ps = chain_s.p_125_all_s
    pb = chain_s.p_125_all_b
    if(ps!=0 and pb!=0):
        ps_ = -ROOT.TMath.Log10(ps)
        pb_ = -ROOT.TMath.Log10(pb)
        
        hs_.Fill(ps_/ (ps_ + pb_) )
                 

print chain_b.GetEntries()
nEntries_b = chain_b.GetEntries()

for i in xrange(nEntries_b):
    chain_b.GetEntry(i)
    ps = chain_b.p_125_all_s
    pb = chain_b.p_125_all_b
    if(ps!=0 and pb!=0):
        ps_ = -ROOT.TMath.Log10(ps)
        pb_ = -ROOT.TMath.Log10(pb)
        hb_.Fill(ps_/ (ps_ + pb_) )


#print hb_.Integral()
#hb_.Scale(1/hb_.Integral())
#hs_.Scale(1/hs_.Integral())
        
hs_.SetTitle("Probability distribution under tHj and ttj hypotheses for bakground events")
hs_.Draw("HIST")
hb_.Draw("HISTSAME")
ms = hs_.GetMaximum()
mb = hb_.GetMaximum()
m = max([ms, mb])
hs_.SetMaximum(m*1.1)


leg.AddEntry(hs_, "tHj events", "f")
leg.AddEntry(hb_, "ttj events", "f")
leg.Draw("SAME")

c.Update()
c.Print("mem_LD.eps")
c.Print("mem_LD.png")    



#ps = np.array([0.])
#tree.SetBranchAddress("p_125_all_s", ps)
#print ps
