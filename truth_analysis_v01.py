###############################################
# Truth analysis 
##############################################

import ROOT

##############################################
# Define the input and output root files
##############################################
f_inputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/mg5_ttbar_jet_001.root", "read")
t = f_inputRoot.Get("ProMC")
type(t)
f_outputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_06_June_15_truth_analysis_v01.root", "recreate")

##############################################
# Define constants
##############################################
m_W = 80.   # GeV 
m_t = 173.  # GeV
nentries = t.GetEntries() # 10000

##############################################
# Define the plots to produce
##############################################
# Plots from input ROOT TFile Tree                 

h_pdg_id = ROOT.TH1F("h_pdg_id", "pdg_id; Number of entries; pgd_id", 2, -4*10**6, 4*10**6)
h_mother1 = ROOT.TH1F("h_mother1", "mother1; Number of entries; mother1", 50, 0, 1200)
h_mother2 = ROOT.TH1F("h_mother2", "mother2; Number of entries; mother2", 50, 0, 1200)

##############################################
# Main analysis - loop over all events
##############################################
for i in range(1000):
    if (( i % 1000 == 0 )): 
       print(":: Processing entry ", i, " = ", i*1.0/nentries*100.0, "%.")    
    if t.LoadTree(i) < 0:
       print("**could not load tree for entry #%s") % i
       break
    nb = t.GetEntry(i) 
    if nb <= 0:
       # no data
       continue
    print('event', i) 
    h_pdg_id.Fill(t.pdg_id) 
    h_mother1.Fill(t.mother1) 
    h_mother2.Fill(t.mother2)  

##############################################
# Draw all histograms and save them.
##############################################
c = ROOT.TCanvas()

h_pdg_id.Draw("E") # put error bars 
c.SaveAs("h_pdg_id.pdf")
h_mother1.Draw("E")
c.SaveAs("h_mother1.pdf")
h_mother2.Draw("E")
c.SaveAs("h_mother2.pdf")

# save to ROOT output files
h_pdg_id.Write()
h_mother1.Write()
h_mother2.Write()

f_outputRoot.Close()