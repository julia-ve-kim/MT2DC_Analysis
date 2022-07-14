#############################################################
# mt2dc final plotting make script
#############################################################

import ROOT
import numpy as np 

##############################################
# Define the input and output root files. 
##############################################
f_inputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_05_May_10_mt2dc_analysis_v01.root", "read")
t = f_inputRoot.Get("results")
type(t)

outDir = "/Users/juliakim/Documents/styledPlotsOutputs/"  
f_outputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_07_July_7_mt2dc_makeFinalPlots.root", "recreate")

##############################################
# Define constants.  
##############################################
m_W = 80 # GeV 
m_t = 173 # GeV
nentries = t.GetEntries() 
num_alpha = len(np.linspace(0, 1, 20))

##############################################
# Produce & fill histograms.  
##############################################
ROOT.gStyle.SetTitleFontSize(0.05)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0)

h2_muT2dc = ROOT.TH2F("h2_muT2dc", "muT2dc; alpha; muT2dc", num_alpha, 0, 1, 30, 0, 3) 
h2_muT2dc_UC = ROOT.TH2F("h2_muT2dc_UC", "muT2dc_UC; alpha; muT2dc_UC", num_alpha, 0, 1, 30, 0, 3) 

h2_muT2prime_W = ROOT.TH2F("h2_muT2prime_W","muT2prime_W; alpha; muT2prime_W", num_alpha, 0, 1, 30, 0, 3)
h2_muT2prime_W_UC = ROOT.TH2F("h2_muT2prime_W_UC","muT2prime_W_UC; alpha; muT2prime_W_UC", num_alpha, 0, 1, 30, 0, 3) 

h2_muT2prime_t = ROOT.TH2F("h2_muT2prime_t","muT2prime_t; alpha; muT2prime_t", num_alpha, 0, 1, 30, 0, 3)
h2_muT2prime_t_UC = ROOT.TH2F("h2_muT2prime_t_UC","muT2prime_t; alpha; muT2prime_t", num_alpha, 0, 1, 30, 0, 3)

h2_mT2dc_diff = ROOT.TH2F("h2_mT2dc_diff","mT2dc_diff; alpha; mT2dc_diff [GeV]", num_alpha, 0, 1, 10, 0, 3*10**4) 
h2_mT2dc_diff_UC = ROOT.TH2F("h2_mT2dc_diff_UC","mT2dc_diff_UC; alpha; mT2dc_diff_UC [GeV]", num_alpha, 0, 1, 10, 0, 3*10**4) 

h2_sub_pT_sideA = ROOT.TH2F("h2_sub_pT_sideA","sub_pT_sideA; alpha; sub_pT_sideA [GeV]", num_alpha, 0, 1, 10, 0, 2*10**7) 
h2_sub_pT_sideA_UC = ROOT.TH2F("h2_sub_pT_sideA_UC","sub_pT_sideA_UC; alpha; sub_pT_sideA_UC [GeV]", num_alpha, 0, 1, 10, 0, 2*10**7) 

h2_sub_pT_sideB = ROOT.TH2F("h2_sub_pT_sideB","sub_pT_sideB; alpha; sub_pT_sideB [GeV]", num_alpha, 0, 1, 10, 0, 2*10**7) 
h2_sub_pT_sideB_UC = ROOT.TH2F("h2_sub_pT_sideB_UC","sub_pT_sideB_UC; alpha; sub_pT_sideB_UC [GeV]", num_alpha, 0, 1, 10, 0, 2*10**7) 

h2_sub_pT_min_over_met = ROOT.TH2F("h2_sub_pT_min_over_met","min(sub_pT_sideA, sub_pT_sideB)/met; alpha; pT_min_over_met [GeV]", num_alpha, 0, 1, 10, 0, 100) 
h2_sub_pT_min_over_met_UC = ROOT.TH2F("h2_sub_pT_min_over_met_UC","min(sub_pT_sideA_UC, sub_pT_sideB_UC)/met; alpha; pT_min_over_met_UC [GeV]", num_alpha, 0, 1, 10, 0, 100) 

#max_W = []
#max_t = []
#max_muT2dc = []
#max_sub_pT_sideA = [] 
#max_mT2dc_diff = []
#max_sub_pT_sideB = []
#max_min_over_met = []

for i in range(nentries): 
    if (( i % 1000 == 0 )): 
       print(":: Processing entry ", i, " = ")    
    if t.LoadTree(i) < 0:
       print("**could not load tree for entry #", i) 
       break
    nb = t.GetEntry(i) 
    if nb <= 0:
       # no data
       continue
    
    if t.constraint_pT_cut == 20 and t.success==True:
        h2_muT2dc.Fill(t.alpha, t.mT2dc/(t.alpha*m_W + (1-t.alpha)*m_t)) 
        h2_muT2prime_W.Fill(t.alpha, t.mT2prime_W/m_W) 
        h2_muT2prime_t.Fill(t.alpha, t.mT2prime_t/m_t) 
        h2_mT2dc_diff.Fill(t.alpha, t.mT2dc_diff)
        h2_sub_pT_sideA.Fill(t.alpha, t.sub_pT_sideA)
        h2_sub_pT_sideB.Fill(t.alpha, t.sub_pT_sideB)
        h2_sub_pT_min_over_met.Fill(t.alpha, t.sub_pT_min_over_met) 
    elif t.constraint_pT_cut == 0 and t.success==True:
        h2_muT2dc_UC.Fill(t.alpha, t.mT2dc/(t.alpha*m_W + (1-t.alpha)*m_t)) 
        h2_muT2prime_W_UC.Fill(t.alpha, t.mT2prime_W/m_W) 
        h2_muT2prime_t_UC.Fill(t.alpha, t.mT2prime_t/m_t) 
        h2_mT2dc_diff_UC.Fill(t.alpha, t.mT2dc_diff)
        h2_sub_pT_sideA_UC.Fill(t.alpha, t.sub_pT_sideA)
        h2_sub_pT_sideB_UC.Fill(t.alpha, t.sub_pT_sideB)
        h2_sub_pT_min_over_met_UC.Fill(t.alpha, t.sub_pT_min_over_met) 

    
##############################################
# Draw all histograms & save as PDFs.
##############################################
c = ROOT.TCanvas()

h2_muT2dc.Draw("COLZ") 
c.SaveAs(outDir+"h2_muT2dc.pdf") 
h2_muT2dc_UC.Draw("COLZ") 
c.SaveAs(outDir+"h2_muT2dc_UC.pdf") 

h2_muT2prime_W.Draw("COLZ") 
c.SaveAs(outDir+"h2_muT2prime_W.pdf")
h2_muT2prime_W_UC.Draw("COLZ") 
c.SaveAs(outDir+"h2_muT2prime_W_UC.pdf") 

h2_muT2prime_t.Draw("COLZ")
c.SaveAs(outDir+"h2_muT2prime_t.pdf") 
h2_muT2prime_t_UC.Draw("COLZ")
c.SaveAs(outDir+"h2_muT2prime_t_UC.pdf") 

h2_mT2dc_diff.Draw("COLZ")
c.SaveAs(outDir+"h2_mT2dc_diff.pdf") 
h2_mT2dc_diff_UC.Draw("COLZ")
c.SaveAs(outDir+"h2_mT2dc_diff_UC.pdf") 

h2_sub_pT_sideA.Draw("COLZ")
c.SaveAs(outDir+"h2_sub_pT_sideA.pdf") 
h2_sub_pT_sideA_UC.Draw("COLZ")
c.SaveAs(outDir+"h2_sub_pT_sideA_UC.pdf") 

h2_sub_pT_sideB.Draw("COLZ")
c.SaveAs(outDir+"h2_sub_pT_sideB.pdf") 
h2_sub_pT_sideB_UC.Draw("COLZ")
c.SaveAs(outDir+"h2_sub_pT_sideB_UC.pdf") 

h2_sub_pT_min_over_met.Draw("COLZ")
c.SaveAs(outDir+"h2_sub_pT_min_over_met.pdf") 
h2_sub_pT_min_over_met_UC.Draw("COLZ")
c.SaveAs(outDir+"h2_sub_pT_min_over_met_UC.pdf") 

##############################################
# Write histograms to output file.
##############################################
h2_muT2dc.Write()
h2_muT2dc_UC.Write()

h2_muT2prime_W.Write()
h2_muT2prime_W_UC.Write()

h2_muT2prime_t.Write() 
h2_muT2prime_t_UC.Write() 

h2_mT2dc_diff.Write()
h2_mT2dc_diff_UC.Write()

h2_sub_pT_sideA.Write()
h2_sub_pT_sideA_UC.Write()

h2_sub_pT_sideB.Write()
h2_sub_pT_sideB_UC.Write()

h2_sub_pT_min_over_met.Write() 
h2_sub_pT_min_over_met_UC.Write() 

f_outputRoot.Close()