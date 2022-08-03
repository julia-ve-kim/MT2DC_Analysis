###############################################
# Truth analysis 
##############################################

import ROOT
import numpy as np 
import mt2dc_v02_prime as DC 

##############################################
# Define the input and output root files
##############################################
f_inputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_06_June_14_truthSkim__mg5_ttbar_jet_001-059_v01.root", "read")
t = f_inputRoot.Get("variables")
type(t)
f_outputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_06_June_15_truth_analysis_v01.root", "recreate")
outDir = "/Users/juliakim/Documents/truthAnalysisPlots/"  

##############################################
# Define constants
##############################################
m_W = 80.   # GeV 
m_t = 173.  # GeV
nentries = t.GetEntries() 

##############################################
# Define the plots to produce
##############################################
# Plots from input ROOT TFile Tree                 

h_W_A_Px = ROOT.TH1F("h_W_A_Px", "W_A_Px; Number of entries / 2 GeV; W_A_Px [GeV]", 300, -600, 600)
h_W_A_Py = ROOT.TH1F("h_W_A_Py", "W_A_Py; Number of entries / 2 GeV; W_A_Py [GeV]", 300, -600, 600)
h_W_A_daughter1 = ROOT.TH1F("h_W_A_daughter1", "W_A_daughter1; Number of entries / 2 units; W_A_daughter1", 40, 40, 120)
h_W_A_daughter2 = ROOT.TH1F("h_W_A_daughter2", "W_A_daughter2; Number of entries / 2 units; W_A_daughter2", 40, 40, 120)
h_W_A_mother1 = ROOT.TH1F("h_W_A_mother1", "W_A_mother1; Number of entries / 2 units; W_A_mother1", 50, 20, 120)
h_W_A_mother2 = ROOT.TH1F("h_W_A_mother2", "W_A_mother2; Number of entries / 2 units; W_A_mother2", 50, 20, 120)
h_W_A_pdg_id = ROOT.TH1F("h_W_A_pdg_id", "W_A_pdg_id; Number of entries / 1 unit; W_A_pdg_id", 4, -26, -22)
h_W_A_M = ROOT.TH1F("h_W_A_M", "W_A_M; Number of entries / 1 unit; W_A_M", 150, 0, 150)

h_ell_A_Px = ROOT.TH1F("h_ell_A_Px", "ell_A_Px; Number of entries / 2 GeV; ell_A_Px [GeV]", 300, -300, 300)
h_ell_A_Py = ROOT.TH1F("h_ell_A_Py", "ell_A_Py; Number of entries / 2 GeV; ell_A_Py [GeV]", 300, -300, 300)
h_ell_A_PT = ROOT.TH1F("h_ell_A_PT", "ell_A_PT; Number of entries / 1 GeV; ell_A_PT [GeV]", 300, 0, 300)
h_ell_A_daughter1 = ROOT.TH1F("h_ell_A_daughter1", "ell_A_daughter1; Number of entries / 2 units; ell_A_daughter1", 40, 40, 120)
h_ell_A_daughter2 = ROOT.TH1F("h_ell_A_daughter2", "ell_A_daughter2; Number of entries / 2 units; ell_A_daughter2", 40, 40, 120)
h_ell_A_mother1 = ROOT.TH1F("h_ell_A_mother1", "ell_A_mother1; Number of entries / 2 units; ell_A_mother1", 50, 20, 120)
h_ell_A_mother2 = ROOT.TH1F("h_ell_A_mother2", "ell_A_mother2; Number of entries / 2 units; ell_A_mother2", 50, 20, 120)
h_ell_A_pdg_id = ROOT.TH1F("h_ell_A_pdg_id", "ell_A_pdg_id; Number of entries / 1 unit; ell_A_pdg_id", 4, 10, 14)

h_ell_B_Px = ROOT.TH1F("h_ell_B_Px", "ell_B_Px; Number of entries / 2 GeV; ell_B_Px [GeV]", 300, -300, 300)
h_ell_B_Py = ROOT.TH1F("h_ell_B_Py", "ell_B_Py; Number of entries / 2 GeV; ell_B_Py [GeV]", 300, -300, 300)
h_ell_B_PT = ROOT.TH1F("h_ell_B_PT", "ell_B_PT; Number of entries / 1 GeV; ell_B_PT [GeV]", 300, 0, 300)
h_ell_B_daughter1 = ROOT.TH1F("h_ell_B_daughter1", "ell_B_daughter1; Number of entries / 2 units; ell_B_daughter1", 40, 30, 110)
h_ell_B_daughter2 = ROOT.TH1F("h_ell_B_daughter2", "ell_B_daughter2; Number of entries / 2 units; ell_B_daughter2", 40, 30, 110)
h_ell_B_mother1 = ROOT.TH1F("h_ell_B_mother1", "ell_B_mother1; Number of entries / 2 units; ell_B_mother1", 40, 20, 100)
h_ell_B_mother2 = ROOT.TH1F("h_ell_B_mother2", "ell_B_mother2; Number of entries / 2 units; ell_B_mother2", 40, 20, 100)
h_ell_B_pdg_id = ROOT.TH1F("h_ell_B_pdg_id", "ell_B_pdg_id; Number of entries / 1 unit; ell_B_pdg_id", 4, -14, -10)
h_W_B_M = ROOT.TH1F("h_W_B_M", "W_B_M; Number of entries / 1 unit; W_B_M", 150, 0, 150)

h_W_A_B_M_max = ROOT.TH1F("h_W_A_B_M_max", "max(W_A_M, W_B_M); Number of entries / 1 unit; max(W_A_M, W_B_M)", 150, 0, 150)

h_nu_A_Px = ROOT.TH1F("h_nu_A_Px", "nu_A_Px; Number of entries / 2 GeV; nu_A_Px [GeV]", 300, -300, 300)
h_nu_A_Py = ROOT.TH1F("h_nu_A_Py", "nu_A_Py; Number of entries / 2 GeV; nu_A_Py [GeV]", 300, -300, 300)
h_nu_A_PT = ROOT.TH1F("h_nu_A_PT", "nu_A_PT; Number of entries / 1 GeV; nu_A_PT [GeV]", 300, 0, 300)
h_nu_A_daughter1 = ROOT.TH1F("h_nu_A_daughter1", "nu_A_daughter1; Number of entries /1 unit; nu_A_daughter1", 1, 0, 1)
h_nu_A_daughter2 = ROOT.TH1F("h_nu_A_daughter2", "nu_A_daughter2; Number of entries / 1 units; nu_A_daughter2", 1, 0, 1)
h_nu_A_mother1 = ROOT.TH1F("h_nu_A_mother1", "nu_A_mother1; Number of entries / 2 units; nu_A_mother1", 40, 40, 120)
h_nu_A_mother2 = ROOT.TH1F("h_nu_A_mother2", "nu_A_mother2; Number of entries / 2 units; nu_A_mother2", 40, 40, 120)
h_nu_A_pdg_id = ROOT.TH1F("h_nu_A_pdg_id", "nu_A_pdg_id; Number of entries / 1 unit; nu_A_pdg_id", 4, -14, -10)

h_nu_B_Px = ROOT.TH1F("h_nu_B_Px", "nu_B_Px; Number of entries / 2 GeV; nu_B_Px [GeV]", 300, -300, 300)
h_nu_B_Py = ROOT.TH1F("h_nu_B_Py", "nu_B_Py; Number of entries / 2 GeV; nu_B_Py [GeV]", 300, -300, 300)
h_nu_B_PT = ROOT.TH1F("h_nu_B_PT", "nu_B_PT; Number of entries / 1 GeV; nu_B_PT [GeV]", 300, 0, 300)

h_nu_min_PT = ROOT.TH1F("h_nu_truth_min_PT", "min(#mbox{pT_{A}}, #mbox{pT_{B}}) [truth]; min(#mbox{pT_{A}}, #mbox{pT_{B}}) [GeV]; Number of entries / 1 GeV", 250, 0, 250)
h_nu_max_PT = ROOT.TH1F("h_nu_truth_min_PT", "min(#mbox{pT_{A}}, #mbox{pT_{B}}) [truth]; min(#mbox{pT_{A}}, #mbox{pT_{B}}) [GeV]; Number of entries / 1 GeV", 300, 0, 300)


# for overlay 
h_max_mT_true = ROOT.TH1F("h_max_mT_true",  "max(#mbox{mT_{1}^{true}}, #mbox{mT_{2}^{true}}); max(#mbox{mT_{1}^{true}}, #mbox{mT_{2}^{true}}) [GeV]; Number of entries / 1 GeV", 150, 0, 150)


##############################################
# Main analysis - loop over all events
##############################################
for i in range(nentries):
    if (( i % 1000 == 0 )): 
       print(":: Processing entry ", i, " = ", i*1.0/nentries*100.0, "%.")    
    if t.LoadTree(i) < 0:
       print("**could not load tree for entry #%s") % i
       break
    nb = t.GetEntry(i) 
    if nb <= 0:
       # no data
       continue
    h_W_A_Px.Fill(t.W_A_Px)
    h_W_A_Py.Fill(t.W_A_Py) 
    h_W_A_daughter1.Fill(t.W_A_daughter1)  
    h_W_A_daughter2.Fill(t.W_A_daughter2) 
    h_W_A_mother1.Fill(t.W_A_mother1) 
    h_W_A_mother2.Fill(t.W_A_mother2)
    h_W_A_pdg_id.Fill(t.W_A_pdg_id) 
    h_W_A_M.Fill(t.W_A_M) 
        
    h_ell_A_Px.Fill(t.ell_A_Px)
    h_ell_A_Py.Fill(t.ell_A_Py)
    h_ell_A_daughter1.Fill(t.ell_A_daughter1) 
    h_ell_A_daughter2.Fill(t.ell_A_daughter2)
    h_ell_A_mother1.Fill(t.ell_A_mother1)
    h_ell_A_mother2.Fill(t.ell_A_mother2)
    h_ell_A_pdg_id.Fill(t.ell_A_pdg_id)
    if np.sqrt(t.ell_A_Px**2 + t.ell_A_Py**2) > 10: 
        h_ell_A_PT.Fill(np.sqrt(t.ell_A_Px**2 + t.ell_A_Py**2)) 

    h_ell_B_Px.Fill(t.ell_B_Px)
    h_ell_B_Py.Fill(t.ell_B_Py)
    h_ell_B_daughter1.Fill(t.ell_B_daughter1) 
    h_ell_B_daughter2.Fill(t.ell_B_daughter2)
    h_ell_B_mother1.Fill(t.ell_B_mother1)
    h_ell_B_mother2.Fill(t.ell_B_mother2)
    h_ell_B_pdg_id.Fill(t.ell_B_pdg_id)
    h_W_B_M.Fill(t.W_B_M)
    if np.sqrt(t.ell_B_Px**2 + t.ell_B_Py**2) > 10: 
        h_ell_B_PT.Fill(np.sqrt(t.ell_B_Px**2 + t.ell_B_Py**2)) 
        
    h_W_A_B_M_max.Fill(max(t.W_A_M, t.W_B_M))
    
    h_nu_A_Px.Fill(t.nu_A_Px)
    h_nu_A_Py.Fill(t.nu_A_Py)
    h_nu_A_PT.Fill(np.sqrt(t.nu_A_Px**2 + t.nu_A_Py**2)) 
    h_nu_A_daughter1.Fill(t.nu_A_daughter1) 
    h_nu_A_daughter2.Fill(t.nu_A_daughter2)
    h_nu_A_mother1.Fill(t.nu_A_mother1)
    h_nu_A_mother2.Fill(t.nu_A_mother2)
    h_nu_A_pdg_id.Fill(t.nu_A_pdg_id)

    h_nu_B_Px.Fill(t.nu_B_Px)
    h_nu_B_Py.Fill(t.nu_B_Py)
    h_nu_B_PT.Fill(np.sqrt(t.nu_B_Px**2 + t.nu_B_Py**2)) 
    
    h_nu_min_PT.Fill(min(np.sqrt(t.nu_B_Px**2 + t.nu_B_Py**2), np.sqrt(t.ell_B_Px**2 + t.ell_B_Py**2))) 
    h_nu_max_PT.Fill(max(np.sqrt(t.nu_B_Px**2 + t.nu_B_Py**2), np.sqrt(t.ell_B_Px**2 + t.ell_B_Py**2)))
    
    # calculate mT1{true}, mT2{true} 
    # mT1 {true}
    
    h_nu_A_array = np.array([t.nu_A_Px, t.nu_A_Py, t.nu_A_Pz, t.nu_A_E]) # invis_sideA_array
    h_ell_A_array = np.array([t.ell_A_Px, t.ell_A_Py, t.ell_A_Pz, t.ell_A_E]) # vis_sideA_array 
    h_nu_B_array = np.array([t.nu_B_Px, t.nu_B_Py, t.nu_A_Pz, t.nu_B_E]) # invis_sideB_array 
    h_ell_B_array = np.array([t.ell_B_Px, t.ell_B_Py, t.ell_B_Pz, t.ell_B_E]) # vis_sideB_array 
    
    mTA_true = DC.mT_arrayCalc(h_ell_A_array, h_nu_A_array) 
    mTB_true = DC.mT_arrayCalc(h_ell_B_array, h_nu_B_array) 
    
    h_max_mT_true.Fill(max(mTA_true, mTB_true))
    

##############################################
# Draw all histograms and save them.
##############################################
c = ROOT.TCanvas()

h_W_A_Px.Draw("E") # put error bars 
c.SaveAs(outDir + "h_W_A_Px.pdf")
h_W_A_Py.Draw("E")
c.SaveAs(outDir + "h_W_A_Py.pdf")
h_W_A_daughter1.Draw("E")
c.SaveAs(outDir + "h_W_A_daughter1.pdf")
h_W_A_daughter2.Draw("E")
c.SaveAs(outDir + "h_W_A_daughter2.pdf")
h_W_A_mother1.Draw("E")
c.SaveAs(outDir + "h_W_A_mother1.pdf")
h_W_A_mother2.Draw("E")
c.SaveAs(outDir + "h_W_A_mother2.pdf")
h_W_A_pdg_id.Draw("E")
c.SaveAs(outDir + "h_W_A_pdg_id.pdf")
h_W_A_M.Draw("E")
c.SaveAs(outDir + "h_W_A_M.pdf") 

h_ell_A_Px.Draw("E")
c.SaveAs(outDir + "h_ell_A_Px.pdf")
h_ell_A_Py.Draw("E")
c.SaveAs(outDir + "h_ell_A_Py.pdf")
h_ell_A_PT.Draw("E")
c.SaveAs(outDir + "h_ell_A_PT.pdf") 
h_ell_A_daughter1.Draw("E")
c.SaveAs(outDir + "h_ell_A_daughter1.pdf")
h_ell_A_daughter2.Draw("E")
c.SaveAs(outDir + "h_ell_A_daughter2.pdf")
h_ell_A_mother1.Draw("E")
c.SaveAs(outDir + "h_ell_A_mother1.pdf")
h_ell_A_mother2.Draw("E")
c.SaveAs(outDir + "h_ell_A_mother2.pdf")
h_ell_A_pdg_id.Draw("E")
c.SaveAs(outDir + "h_ell_A_pdg_id.pdf")

h_ell_B_Px.Draw("E")
c.SaveAs(outDir + "h_ell_B_Px.pdf")
h_ell_B_Py.Draw("E")
c.SaveAs(outDir + "h_ell_B_Py.pdf")
h_ell_B_PT.Draw("E")
c.SaveAs(outDir + "h_ell_B_PT.pdf") 
h_ell_B_daughter1.Draw("E")
c.SaveAs(outDir + "h_ell_B_daughter1.pdf")
h_ell_B_daughter2.Draw("E")
c.SaveAs(outDir + "h_ell_B_daughter2.pdf")
h_ell_B_mother1.Draw("E")
c.SaveAs(outDir + "h_ell_B_mother1.pdf")
h_ell_B_mother2.Draw("E")
c.SaveAs(outDir + "h_ell_B_mother2.pdf")
h_ell_B_pdg_id.Draw("E")
c.SaveAs(outDir + "h_ell_B_pdg_id.pdf")
h_W_B_M.Draw("E")
c.SaveAs(outDir + "h_W_B_M.pdf") 

h_W_A_B_M_max.Draw("E")
c.SaveAs(outDir + "h_W_A_B_M_max.pdf") 

h_nu_A_Px.Draw("E")
c.SaveAs(outDir + "h_nu_A_Px.pdf")
h_nu_A_Py.Draw("E")
c.SaveAs(outDir + "h_nu_A_Py.pdf")
h_nu_A_PT.Draw("E")
c.SaveAs(outDir + "h_nu_A_PT.pdf") 
h_nu_A_daughter1.Draw("E")
c.SaveAs(outDir + "h_nu_A_daughter1.pdf")
h_nu_A_daughter2.Draw("E")
c.SaveAs(outDir + "h_nu_A_daughter2.pdf")
h_nu_A_mother1.Draw("E")
c.SaveAs(outDir + "h_nu_A_mother1.pdf")
h_nu_A_mother2.Draw("E")
c.SaveAs(outDir + "h_nu_A_mother2.pdf")
h_nu_A_pdg_id.Draw("E")
c.SaveAs(outDir + "h_nu_A_pdg_id.pdf")

h_nu_B_Px.Draw("E")
c.SaveAs(outDir + "h_nu_B_Px.pdf")
h_nu_B_Py.Draw("E")
c.SaveAs(outDir + "h_nu_B_Py.pdf")
h_nu_B_PT.Draw("E")
c.SaveAs(outDir + "h_nu_B_PT.pdf") 

h_nu_min_PT.Draw("E")
c.SaveAs(outDir + "h_nu_min_PT.pdf") 

h_nu_max_PT.Draw("E")
c.SaveAs(outDir + "h_nu_max_PT.pdf") 

h_max_mT_true.Draw("E")
c.SaveAs(outDir + "h_max_mT_true.pdf") 

# save to ROOT output files
h_W_A_Px.Write()
h_W_A_Py.Write()
h_W_A_daughter1.Write()
h_W_A_daughter2.Write() 
h_W_A_mother1.Write() 
h_W_A_mother2.Write() 
h_W_A_pdg_id.Write() 

h_ell_A_Px.Write()
h_ell_A_Py.Write()
h_ell_A_PT.Write() 
h_ell_A_daughter1.Write()
h_ell_A_daughter2.Write()
h_ell_A_mother1.Write()
h_ell_A_mother2.Write()
h_ell_A_pdg_id.Write()

h_ell_B_Px.Write()
h_ell_B_Py.Write()
h_ell_B_PT.Write() 
h_ell_B_daughter1.Write()
h_ell_B_daughter2.Write()
h_ell_B_mother1.Write()
h_ell_B_mother2.Write()
h_ell_B_pdg_id.Write()

h_W_A_B_M_max.Write()

h_nu_A_Px.Write() 
h_nu_A_Py.Write() 
h_nu_A_PT.Write() 
h_nu_A_daughter1.Write() 
h_nu_A_daughter2.Write() 
h_nu_A_mother1.Write() 
h_nu_A_mother2.Write() 
h_nu_A_pdg_id.Write() 

h_nu_B_Px.Write() 
h_nu_B_Py.Write() 
h_nu_B_PT.Write() 

h_nu_min_PT.Write() 
h_nu_max_PT.Write() 

h_max_mT_true.Write() 

f_outputRoot.Close()