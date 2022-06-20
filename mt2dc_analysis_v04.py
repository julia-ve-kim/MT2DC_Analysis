###############################################
# Mt2dc analysis
# Purpose: To analysis the results of the m2dc
#    variable calculation to determine its properties,
#    it benefits, and its drawbacks.  Comparisons are
#    made with the original mt2 variable.  Final stylized
#    plots are not produced here but by mt2dc_makeFinalPlots.py.
#    This code executes the mt2dc calculation, performed 
#    by mt2dc.py.
##############################################

import ROOT
import mt2dc_v02_prime as DC 
import math
import numpy as np 
import scipy.optimize as so 

##############################################
# Define the input and output root files
##############################################
f_inputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_03_March_07_skim_mg5_ttbar_jet_merged_001-716_ntuple_2l2b_v01.root", "read")
t = f_inputRoot.Get("variables")
type(t)
f_outputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_05_May_10_mt2dc_analysis_v01.root", "recreate")

##############################################
# Define the plots to produce
##############################################
# Plots from input ROOT TFile Tree                 

h_ell1_pt = ROOT.TH1F("h_ell1_pt", "Pt of highest pt light lepton; Leading light lepton pt [GeV]; Number of entries / 2 GeV", 100, 0, 200)
h_ell1_E = ROOT.TH1F("h_ell1_E", "E of highest pt light lepton; Leading light lepton Energy [GeV]; Number of entries / 2.5 GeV", 100, 0, 250)

h_ell2_pt = ROOT.TH1F("h_ell2_pt", "Pt of lowest pt light lepton; Second light lepton pt [GeV]; Number of entries / 2GeV", 100, 0, 200)
h_ell2_E = ROOT.TH1F("h_ell2_pt", "E of lowest pt light lepton; Second light lepton pt [GeV]; Number of entries / 2.5eV", 100, 0, 250)

h_bjet1_E = ROOT.TH1F("h_bjet1_E", "E of highest pt b-tagged jet; Leading b-tagged jet Energy [GeV];Number of entries / 5 GeV", 100, 0, 500)
h_bjet2_E = ROOT.TH1F("h_bjet_2_E", "E of lowest pt b-tagged jet; Leading b-tagged jet Energy [GeV]; Number of entries / 5 GeV", 100, 0, 500)

h_mT2_W = ROOT.TH1F("h_mT2_W", "mt2(ell1,ell2) = mt2(W); mt2(W) [GeV]; Number of entries / 1 GeV", 200, 0, 200)
h_mT2_t_11_22 = ROOT.TH1F("h_mT2_t_11_22", "mt2|11,22(t) = mt2(b1 ell1,b2 ell2); mt2(t|11,22) [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_mT2_t_12_21 = ROOT.TH1F("h_mT2_t_12_21", "mt2|12,21(t) = mt2(b1 ell2,b2 ell1); mt2(t|12,21) [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_mT2_t_min   = ROOT.TH1F("h_mT2_t_min", "min( mt2(t|11,22), mt2(t|12,21)); mt2(t)_min [GeV]; Number of entries / 1 GeV",     300, 0, 300)

h_mT_ell1__forMt2Overlay = ROOT.TH1F("h_mT_ell1_forMt2Overlay", "mt(ell1, EtMiss); mt(ell1) [GeV]; Number of entries / 1 GeV",   300, 0, 300)
h_mT_ell1 = ROOT.TH1F("h_mT_ell1", "mt(ell1, EtMiss); mt(ell1) [GeV]; Number of entries / 5 GeV", 200, 0, 1000)
h_nu_ell1_pt = ROOT.TH1F("h_nu_ell1_pt", "True pt of neutrino ass't with l1; Energy [GeV]; Number of entries / 1GeV", 300, 0, 300)
h_nu_ell1_E = ROOT.TH1F("h_nu_ell1_E", "True E of neutrino ass't with l1; Energy [GeV]; Number of entries / 1GeV", 300, 0, 300)

h_nu_ell2_pt = ROOT.TH1F("h_nu_ell2_pt", "True pt of neutrino ass't with l2; Energy [GeV]; Number of entries / 1GeV", 300, 0, 300)
h_nu_ell2_E = ROOT.TH1F("h_nu_ell2_E", "True E of neutrino ass't with l2; Energy [GeV]; Number of entries / 1GeV", 300, 0, 300)

h_EtMiss = ROOT.TH1F("h_EtMiss", "Missing transverse E; Energy [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_EtMiss_phi = ROOT.TH1F("h_EtMiss_phi", "Azimuthal direction of missing transverse E; Azimuthal angle [rad]; Number of entries / 0.02 rad", 400, -4, 4)  

h_eventNumber = ROOT.TH1F("h_event_num", "Event Number; Event number; Number of enthries / 30E8", 10, 10**9, 40*10**9)
h_truth_numNu = ROOT.TH1F("h_truth_numNu", "True number of ejected neutrinos; Number; Events / 1 u ", 10, 0, 10) 

# Plots of mt2dc calculation   
# alpha = 0 (UC = unconstrained) 
h_mT2dc_diff_alpha_0  = ROOT.TH1F("h_mT2dc_diff_alpha_0", "mT2dc(alpha = 0) - mt2_t_bjet1ell1_bjet2ell2; Difference [GeV]; Number of entries / 2 GeV", 100, -100, 100)
h_mT2dc_alpha_0 = ROOT.TH1F("h_mT2dc_alpha_0", "mT2dc(alpha = 0); mT2dc [GeV]; Number of entries / 3 GeV", 100, 0, 300)
h_mT2prime_W_alpha_0 = ROOT.TH1F("h_mT2prime_W_alpha_0", "mt2'(W); mt2'(W) [GeV]; Number of entries / 3 GeV", 300, 0, 300)
h_mT2prime_t_alpha_0 = ROOT.TH1F("h_mT2prime_t_alpha_0", "mt2'(t); mt2'(t) [GeV]; Number of entries / 3 GeV", 100, 0, 300)
h_pT_alpha_0 = ROOT.TH1F("h_pT_alpha_0", "pT_alpha_0(t); pT [GeV]; Number of entries / 1 GeV", 300, 0, 300)

h_mT2dc_diff_alpha_0_UC  = ROOT.TH1F("h_mT2dc_diff_alpha_0_UC", "mT2dc(alpha = 0) - mt2_t_bjet1ell1_bjet2ell2 [no constraint]; Difference [GeV]; Number of entries / 2 GeV", 100, -100, 100)
h_mT2dc_alpha_0_UC = ROOT.TH1F("h_mT2dc_alpha_0_UC", "mT2dc(alpha = 0) [no constraint]; mT2dc [GeV]; Number of entries / 3 GeV", 100, 0, 300)
h_mT2prime_W_alpha_0_UC = ROOT.TH1F("h_mT2prime_W_alpha_0_UC", "mt2'(W) [no constraint]; mt2'(W) [GeV]; Number of entries / 3 GeV", 300, 0, 300)
h_mT2prime_t_alpha_0_UC = ROOT.TH1F("h_mT2prime_t_alpha_0_UC", "mt2'(t) [no constraint]; mt2'(t) [GeV]; Number of entries / 3 GeV", 100, 0, 300)
h_pT_sideA_alpha_0_UC = ROOT.TH1F("h_pT_sideA_alpha_0_UC", "pT_sideA_alpha_0 [no constraint]; pT [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_pT_sideB_alpha_0_UC = ROOT.TH1F("h_pT_sideB_alpha_0_UC", "pT_sideB_alpha_0 [no constraint]; pT [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_pT_min_alpha_0_UC = ROOT.TH1F("h_pT_min_alpha_0_UC", "min(pT_sideA_alpha_0, pT_sideB_alpha_0) [no constraint]; pT [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_pT_ratio_alpha_0_UC = ROOT.TH1F("h_pT_ratio_alpha_0_UC", "|pT_sideA|/|met| (alpha = 0) [no constraint]; |pT_sideA|/|met|; Number of entries / 0.05 ", 20, 0, 1)

# alpha = 1 (UC = unconstrained) 
h_mT2dc_diff_alpha_1 = ROOT.TH1F("h_mT2dc_diff_alpha_1", "mT2dc(alpha = 1) - mT2(W); Difference [GeV]; Number of entries / 2 GeV", 100, -100, 100)
h_mT2dc_alpha_1 = ROOT.TH1F("h_mT2dc_alpha_1", "mT2dc(alpha = 1); mT2dc [GeV]; Number of entries / 1 GeV", 200, 0, 200)
h_mT2prime_W_alpha_1 = ROOT.TH1F("h_mT2prime_W_alpha_1", "mt2'(W); mt2'(W) [GeV]; Number of entries / 2 GeV", 100, 0, 200)
h_mT2prime_t_alpha_1 = ROOT.TH1F("h_mT2prime_t_alpha_1", "mt2'(t); mt2'(t) [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_pT_alpha_1 = ROOT.TH1F("h_pT_alpha_1", "pT_alpha_1(t); pT [GeV]; Number of entries / 1 GeV", 300, 0, 300)

h_mT2dc_diff_alpha_1_UC = ROOT.TH1F("h_mT2dc_diff_alpha_1_UC", "mT2dc(alpha = 1) - mT2(W) [no constraint]; Difference [GeV]; Number of entries / 2 GeV", 100, -100, 100)
h_mT2dc_alpha_1_UC = ROOT.TH1F("h_mT2dc_alpha_1_UC", "mT2dc(alpha = 1) [no constraint]; mT2dc [GeV]; Number of entries / 1 GeV", 200, 0, 200)
h_mT2prime_W_alpha_1_UC = ROOT.TH1F("h_mT2prime_W_alpha_1_UC", "mt2'(W) [no constraint]; mt2'(W) [GeV]; Number of entries / 2 GeV", 100, 0, 200)
h_mT2prime_t_alpha_1_UC = ROOT.TH1F("h_mT2prime_t_alpha_1_UC", "mt2'(t) [no constraint]; mt2'(t) [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_pT_sideA_alpha_1_UC = ROOT.TH1F("h_pT_sideA_alpha_1_UC", "pT_sideA_alpha_1 [no constraint]; pT [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_pT_sideB_alpha_1_UC = ROOT.TH1F("h_pT_sideB_alpha_1_UC", "pT_sideB_alpha_1 [no constraint]; pT [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_pT_min_alpha_1_UC = ROOT.TH1F("h_pT_min_alpha_1_UC", "min(pT_sideA_alpha_1, pT_sideB_alpha_1) [no constraint]; pT [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_pT_ratio_alpha_1_UC = ROOT.TH1F("h_pT_ratio_alpha_1_UC", "|pT_sideA|/|met| (alpha = 1) [no constraint]; |pT_sideA|/|met|; Number of entries / 0.05 ", 20, 0, 1)

# 2-D Histogram 
ROOT.gStyle.SetTitleFontSize(0.05)
ROOT.gStyle.SetLabelFont(42, "XYZ")
ROOT.gStyle.SetTitleOffset(1)
ROOT.gStyle.SetTitleOffset(2, "XY")
ROOT.gStyle.SetTitleOffset(1.25, "Z")

h_mT2dc_alpha_0_C_vs_UC = ROOT.TH2F("h_mT2dc_alpha_0_C_vs_UC", "(mT2dc(t) [no constraint], mT2dc(t) [constraint]); mT2dc(t) [no constraint]; mT2dc(t) [constraint]; Number of entries / (4, 4) GeV ", 65, 40, 300, 40, 30, 300) 
h_mT2dc_alpha_1_C_vs_UC = ROOT.TH2F("h_mT2dc_alpha_1_C_vs_UC", "(mT2dc(W) [no constraint], mT2dc(W) [constraint]); mT2dc(W) [no constraint]; mT2dc(W) [constraint]; Number of entries / (4, 4) GeV ", 25, 0, 100, 25, 0, 100) 

h_alpha_1_mT2dcW_vs_minPT = ROOT.TH2F("h_alpha_1_mT2dcW_vs_minPT", "(mT2dc(W) [constraint], min(pT)); mT2dc(W) [constraint] [GeV]; min(pT_sideA_alpha_1, pT_sideB_alpha_1) [GeV]; Number of entries / (1, 1) GeV ", 100, 0, 100, 100, 0, 100) 


##############################################
# Define constants
##############################################
m_W = 80.   # GeV 
m_t = 173.  # GeV
nentries = t.GetEntries() # 60599
alphaListindex = input('Choose 0 (alpha = 0) or 1 (alpha = 1):') 
calcStyle = input('Choose fast or slow:') 
what_index = [] 

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
    
    # REQUIRED IN MINIMISATION: Retrive information from tree & successively fill histograms. 
    # get sideA bjet information (suppose sideA = highest pt jet originating from a b-quark)
    bjet1_sideA_Px = DC.extract_Px(t.bjet1_PT, t.bjet1_Eta, t.bjet1_Phi, t.bjet1_Mass) 
    bjet1_sideA_Py = DC.extract_Py(t.bjet1_PT, t.bjet1_Eta, t.bjet1_Phi, t.bjet1_Mass)
    bjet1_sideA_Pz = DC.extract_Pz(t.bjet1_PT, t.bjet1_Eta, t.bjet1_Phi, t.bjet1_Mass)
    bjet1_sideA_E = DC.extract_E(t.bjet1_PT, t.bjet1_Eta, t.bjet1_Phi, t.bjet1_Mass) 
    bjet1_sideA_array = np.array([bjet1_sideA_Px, bjet1_sideA_Py, bjet1_sideA_Pz, bjet1_sideA_E]) 

    h_bjet1_E.Fill(bjet1_sideA_E)
    
    # get sideA lepton information (suppose sideA = highest p light lepton) 
    ell1_sideA_Px = DC.extract_Px(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0) 
    ell1_sideA_Py = DC.extract_Py(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0)
    ell1_sideA_Pz = DC.extract_Pz(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0)
    ell1_sideA_E = DC.extract_E(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0) 
    ell1_sideA_array = np.array([ell1_sideA_Px, ell1_sideA_Py, ell1_sideA_Pz, ell1_sideA_E]) 
    
    h_ell1_pt.Fill(np.sqrt(ell1_sideA_Px**2 + ell1_sideA_Px**2)) 
    h_ell1_E.Fill(ell1_sideA_E) 
    
    vis_sideA_array = np.array([bjet1_sideA_array, ell1_sideA_array]) 
    
    # get sideB bjet information
    bjet2_sideB_Px = DC.extract_Px(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass) 
    bjet2_sideB_Py = DC.extract_Py(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass)
    bjet2_sideB_Pz = DC.extract_Pz(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass)
    bjet2_sideB_E = DC.extract_E(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass) 
    bjet2_sideB_array = np.array([bjet2_sideB_Px, bjet2_sideB_Py, bjet2_sideB_Pz, bjet2_sideB_E])
    
    h_bjet2_E.Fill(bjet2_sideB_E)
    
    # get sideB lepton information (suppose sideB = second highest p light lepton) 
    ell2_sideB_Px = DC.extract_Px(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0) 
    ell2_sideB_Py = DC.extract_Py(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0)
    ell2_sideB_Pz = DC.extract_Pz(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0)
    ell2_sideB_E = DC.extract_E(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0) 
    ell2_sideB_array = np.array([ell2_sideB_Px, ell2_sideB_Py, ell2_sideB_Pz, ell2_sideB_E]) 
        
    h_ell2_pt.Fill(np.sqrt(ell2_sideB_Px**2 + ell2_sideB_Py**2))
    h_ell2_E.Fill(ell2_sideB_E)
   
    vis_sideB_array = np.array([bjet2_sideB_array, ell2_sideB_array]) 
    
    # get met information 
    met_Px = DC.extract_Px(t.EtMiss, 0, t.EtMiss_phi, 0) 
    met_Py = DC.extract_Py(t.EtMiss, 0, t.EtMiss_phi, 0) 
    met_E = DC.extract_E(t.EtMiss, 0, t.EtMiss_phi, 0) 
    met = np.array([met_Px, met_Py, 0, met_E])
    
    mt_ell1 = DC.mT_arrayCalc(ell1_sideA_array, met)
    h_mT_ell1.Fill(mt_ell1)
    h_mT_ell1__forMt2Overlay.Fill(mt_ell1)
    
    h_EtMiss.Fill(t.EtMiss) 
    h_EtMiss_phi.Fill(t.EtMiss_phi) 
    
    # get transverse mass information
    mt2_W = t.mt2_W_ell1ell2
    mt2_t_11_22 = t.mt2_t_bjet1ell1_bjet2ell2
    mt2_t_12_21 = t.mt2_t_bjet1ell2_bjet2ell1
    mt2_t_min = min(mt2_t_11_22, mt2_t_12_21)

    h_mT2_W.Fill(mt2_W)
    h_mT2_t_11_22.Fill(mt2_t_11_22)
    h_mT2_t_12_21.Fill(mt2_t_12_21)
    h_mT2_t_min.Fill(mt2_t_min)
    
    # NOT REQUIRED IN MINIMISATION: Retrive information from tree & successively fill histograms. 
    # Truth neutrino 1 
    pt  = t.truth_nu_ell1_PT
    eta = t.truth_nu_ell1_Eta
    phi = t.truth_nu_ell1_Phi
    m = 0
    p4_nu_ell1 = ROOT.TLorentzVector() 
    p4_nu_ell1.SetPtEtaPhiM(pt,eta,phi,m)

    h_nu_ell1_pt.Fill(p4_nu_ell1.Pt()) 
    h_nu_ell1_E.Fill(p4_nu_ell1.E()) 
    
    # Truth neutrino 2 
    pt  = t.truth_nu_ell2_PT
    eta = t.truth_nu_ell2_Eta
    phi = t.truth_nu_ell1_Phi
    m = 0
    p4_nu_ell2 = ROOT.TLorentzVector() 
    p4_nu_ell2.SetPtEtaPhiM(pt,eta,phi,m)

    h_nu_ell2_pt.Fill(p4_nu_ell2.Pt()) 
    h_nu_ell2_E.Fill(p4_nu_ell2.E()) 
    
    # True number of neutrinos 
    truth_numNu = t.truth_numNu
    h_truth_numNu.Fill(truth_numNu) 

    # Event number 
    event_num = t.eventNumber
    h_eventNumber.Fill(event_num)
    
    # CALCULATION OF MT2DC:
    alphaList = [0, 1] 
    index = int(alphaListindex) 
    invis_sideA_array_guess = met[:2]/2 
    invis_sideA_array_guess_2 = bjet1_sideA_array[:2] 
    invis_sideA_array_guess_3 = ell1_sideA_array[:2] 
    invis_sideA_array_guess_4 = bjet2_sideB_array[:2] 
    invis_sideA_array_guess_5 = ell2_sideB_array[:2] 
    invis_sideA_array_guess_6 = 1.01*met[:2] 
    invis_sideA_array_guess_7 = 0.01*met[:2] 
    invis_sideA_array_guess_8 = 2*met[:2] 
    invis_sideA_array_guess_9 = ell1_sideA_array[:2] + bjet1_sideA_array[:2] 
    invis_sideA_array_guess_10 = ell2_sideB_array[:2] + bjet2_sideB_array[:2] 
    invis_sideA_array_guess_11 = ell1_sideA_array[:2] + bjet2_sideB_array[:2] 
    invis_sideA_array_guess_12 = ell2_sideB_array[:2] + bjet1_sideA_array[:2] 
    invis_sideA_array_guess_13 = met[:2] + ell1_sideA_array[:2]
    invis_sideA_array_guess_14 = met[:2] + ell2_sideB_array[:2]   
    invis_sideA_array_guess_15 = met[:2] + bjet1_sideA_array[:2]
    invis_sideA_array_guess_16 = met[:2] + bjet2_sideB_array[:2]   
    invis_sideA_array_guess_17 = -met[:2] 
    invis_sideA_array_guess_18 = met[:2] - ell1_sideA_array[:2]
    invis_sideA_array_guess_19 = met[:2] - ell2_sideB_array[:2]
    invis_sideA_array_guess_20 = met[:2] - bjet1_sideA_array[:2]  
    invis_sideA_array_guess_21 = met[:2] - bjet2_sideB_array[:2]
    
    def objective(invis_sideA_2vec): # minimise over a 2-vector array, having components px and py 
        invis_sideA_array = np.array([invis_sideA_2vec[0], invis_sideA_2vec[1], 0, 
                                     np.sqrt(invis_sideA_2vec[0]**2 + invis_sideA_2vec[0]**2)]) 
        
        alpha_term_1 = DC.mT_arrayCalc(vis_sideA_array[-1], invis_sideA_array) # mT(lA, pT_A)
        alpha_term_2 = DC.mT_arrayCalc(vis_sideB_array[-1], met-invis_sideA_array) # mT(TB, pT_B) 
        alpha_term = max(alpha_term_1, alpha_term_2) 
        
        beta_term_1 = DC.mT_arrayCalc(vis_sideA_array[0] + vis_sideA_array[-1], invis_sideA_array) # mT(lATA, pT_A)
        beta_term_2 = DC.mT_arrayCalc(vis_sideB_array[0] + vis_sideB_array[-1], met-invis_sideA_array) # mT(TBbB, pt_B)
        beta_term = max(beta_term_1, beta_term_2) 
    
        return alphaList[index]*alpha_term + (1-alphaList[index])*beta_term 
    
    pT_lower_bound = 10
    def constraint_1(invis_sideA_2vec):
        return invis_sideA_2vec[0]**2 + invis_sideA_2vec[1]**2 - pT_lower_bound**2 
    
    def constraint_2(invis_sideA_2vec):
        invis_sideB_2vec = met[:2] - invis_sideA_2vec 
        return invis_sideB_2vec[0]**2 + invis_sideB_2vec[1]**2 - pT_lower_bound**2 
        
    if calcStyle == 'fast': 
        invis_sideA_array_guesses = [invis_sideA_array_guess, invis_sideA_array_guess_2, invis_sideA_array_guess_3, 
                                     invis_sideA_array_guess_4, invis_sideA_array_guess_5, invis_sideA_array_guess_6, 
                                     invis_sideA_array_guess_7, invis_sideA_array_guess_8, invis_sideA_array_guess_9, 
                                     invis_sideA_array_guess_10, invis_sideA_array_guess_11, invis_sideA_array_guess_12,
                                     invis_sideA_array_guess_13, invis_sideA_array_guess_14, invis_sideA_array_guess_15,
                                     invis_sideA_array_guess_16, invis_sideA_array_guess_17, invis_sideA_array_guess_18, 
                                    invis_sideA_array_guess_19, invis_sideA_array_guess_20, invis_sideA_array_guess_21]
        guess = objective(invis_sideA_array_guess)
        guess_2 = objective(invis_sideA_array_guess_2) 
        guess_3 = objective(invis_sideA_array_guess_3) 
        guess_4 = objective(invis_sideA_array_guess_4)
        guess_5 = objective(invis_sideA_array_guess_5) 
        guess_6 = objective(invis_sideA_array_guess_6)
        guess_7 = objective(invis_sideA_array_guess_7)
        guess_8 = objective(invis_sideA_array_guess_8)
        guess_9 = objective(invis_sideA_array_guess_9)        
        guess_10 = objective(invis_sideA_array_guess_10)
        guess_11 = objective(invis_sideA_array_guess_11)        
        guess_12 = objective(invis_sideA_array_guess_12)           
        guess_13 = objective(invis_sideA_array_guess_13)  
        guess_14 = objective(invis_sideA_array_guess_14)         
        guess_15 = objective(invis_sideA_array_guess_15)  
        guess_16 = objective(invis_sideA_array_guess_16)          
        guess_17 = objective(invis_sideA_array_guess_17)   
        guess_18 = objective(invis_sideA_array_guess_18)   
        guess_19 = objective(invis_sideA_array_guess_19)   
        guess_20 = objective(invis_sideA_array_guess_20)   
        guess_21 = objective(invis_sideA_array_guess_21)
                
        guesses = [guess, guess_2, guess_3, guess_4, guess_5, guess_6, guess_7, guess_8, guess_9, guess_10, guess_11,
                  guess_12, guess_13, guess_14, guess_15, guess_16, guess_17, guess_18, guess_19, guess_20, guess_21] 
        
        cons = [{'type': 'ineq', 'fun': constraint_1},  {'type': 'ineq', 'fun': constraint_2}] 
        
        sol_UC = so.minimize(objective, x0 = invis_sideA_array_guesses[np.argmin(guesses)], method='SLSQP', 
                          options={'maxiter': 2000, 'ftol': 1e-07,'disp': True})  
        sol = so.minimize(objective, x0 = invis_sideA_array_guesses[np.argmin(guesses)], method='SLSQP', 
                          options={'maxiter': 2000, 'ftol': 1e-07,'disp': True}, constraints=cons)  
        
        print('np.argmin(guesses)', np.argmin(guesses))
        what_index.append([np.argmin(guesses)]) 
        
        if index == 0:
            h_mT2dc_diff_alpha_0.Fill(sol.fun - mt2_t_11_22) 
            h_mT2dc_alpha_0.Fill(sol.fun) 
            h_mT2prime_W_alpha_0.Fill(DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol.x)) 
            h_mT2prime_t_alpha_0.Fill(DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol.x))  
            h_pT_alpha_0.Fill(np.linalg.norm(sol.x))

            # unconstrained 
            h_mT2dc_diff_alpha_0_UC.Fill(sol_UC.fun - mt2_t_11_22) 
            h_mT2dc_alpha_0_UC.Fill(sol_UC.fun) 
            h_mT2prime_W_alpha_0_UC.Fill(DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x)) 
            h_mT2prime_t_alpha_0_UC.Fill(DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x))  
            h_pT_sideA_alpha_0_UC.Fill(np.linalg.norm(sol_UC.x)) # magnitude pT A side 
            h_pT_sideB_alpha_0_UC.Fill(np.linalg.norm(met[:2] - sol_UC.x)) # magnitude pT B side 
            h_pT_min_alpha_0_UC.Fill(min(np.linalg.norm(sol_UC.x), np.linalg.norm(met[:2] - sol_UC.x)))
            h_pT_ratio_alpha_0_UC.Fill(np.linalg.norm(sol_UC.x)/np.linalg.norm(met[:2])) 
            
            h_mT2dc_alpha_0_C_vs_UC.Fill(sol_UC.fun, sol.fun)
            
        
        elif index == 1:
            h_mT2dc_diff_alpha_1.Fill(sol.fun - mt2_W)  
            h_mT2dc_alpha_1.Fill(sol.fun) 
            h_mT2prime_W_alpha_1.Fill(DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol.x))  
            h_mT2prime_t_alpha_1.Fill(DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol.x))  
            h_pT_alpha_1.Fill(np.linalg.norm(sol.x))
            
            h_alpha_1_mT2dcW_vs_minPT.Fill(sol.fun, min(np.linalg.norm(sol.x), np.linalg.norm(met[:2] - sol.x))) 

            # unconstrained 
            h_mT2dc_diff_alpha_1_UC.Fill(sol_UC.fun - mt2_W)  
            h_mT2dc_alpha_1_UC.Fill(sol_UC.fun) 
            h_mT2prime_W_alpha_1_UC.Fill(DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x))  
            h_mT2prime_t_alpha_1_UC.Fill(DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x))  
            h_pT_sideA_alpha_1_UC.Fill(np.linalg.norm(sol_UC.x))
            h_pT_sideB_alpha_1_UC.Fill(np.linalg.norm(met[:2] - sol_UC.x))           
            h_pT_min_alpha_1_UC.Fill(min(np.linalg.norm(sol_UC.x), np.linalg.norm(met[:2] - sol_UC.x))) 
            h_pT_ratio_alpha_1_UC.Fill(np.linalg.norm(sol_UC.x)/np.linalg.norm(met[:2])) 
            
            h_mT2dc_alpha_1_C_vs_UC.Fill(sol_UC.fun, sol.fun)
    
    elif calcStyle == 'slow':
        sol_1 = so.minimize(objective, x0 = invis_sideA_array_guess_1, method='SLSQP', 
                            options={'maxiter': 2000, 'ftol': 1e-07,'disp': True}) 
        sol_2 = so.minimize(objective, x0 = invis_sideA_array_guess_2, method='SLSQP', 
                        options={'maxiter': 2000, 'ftol': 1e-07,'disp': True}) 
        sol_3 = so.minimize(objective, x0 = invis_sideA_array_guess_3, method='SLSQP', 
                        options={'maxiter': 2000, 'ftol': 1e-07,'disp': True})     
        sol_4 = so.minimize(objective, x0 = invis_sideA_array_guess_4, method='SLSQP', 
                        options={'maxiter': 2000, 'ftol': 1e-07,'disp': True}) 
        sol_5 = so.minimize(objective, x0 = invis_sideA_array_guess_5, method='SLSQP', 
                        options={'maxiter': 2000, 'ftol': 1e-07,'disp': True})
        sol_6 = so.minimize(objective, x0 = invis_sideA_array_guess_6, method='SLSQP', 
                        options={'maxiter': 2000, 'ftol': 1e-07,'disp': True})  
        sol_7 = so.minimize(objective, x0 = invis_sideA_array_guess_7, method='SLSQP', 
                        options={'maxiter': 2000, 'ftol': 1e-07,'disp': True}) 
        
        sol_fun_array = [sol_1.fun, sol_2.fun, sol_3.fun, sol_4.fun, sol_5.fun, sol_6.fun, sol_7.fun] 
        sol_fun = min(sol_fun_array) 
        sol_x_array = [sol_1.x, sol_2.x, sol_3.x, sol_4.x, sol_5.x, sol_6.x, sol_7.x]
        sol_x = sol_x_array[np.argmin(sol_fun_array)] 
        
        if index == 0:
            h_mT2dc_diff_alpha_0.Fill(sol_fun - mt2_t_11_22) 
            h_mT2dc_alpha_0.Fill(sol_fun) 
            h_mT2prime_t_alpha_0.Fill(DC.get_beta_term(vis_sideA_array, vis_sideB_aray, met, sol_x))  
        elif index == 1:
            h_mT2dc_diff_alpha_1.Fill(sol_fun - mt2_W) 
            h_mT2dc_alpha_1.Fill(sol_fun) 
            h_mT2prime_W_alpha_1.Fill(DC.get_alpha_term(vis_sideA_array, vis_sideB_aray, met, sol_x))  
            
print(np.mean(what_index), np.std(what_index, ddof=1)) 
print(what_index) 
            
##############################################
# Draw all histograms and save them.
##############################################
c = ROOT.TCanvas()

h_ell1_pt.Draw("E") # put error bars 
c.SaveAs("h_ell1_PT.pdf")
h_ell1_E.Draw("E")
c.SaveAs("h_ell1_E.pdf")
h_ell2_pt.Draw("E")
c.SaveAs("h_ell2_PT.pdf")
h_ell2_E.Draw("E")
c.SaveAs("h_ell2_E.pdf")

h_bjet1_E.Draw("E")
c.SaveAs("h_bjet1_E.pdf")
h_bjet2_E.Draw("E")
c.SaveAs("h_bjet2_E.pdf")

h_nu_ell1_E.Draw("E")
c.SaveAs("h_nu_ell1_E.pdf")
h_nu_ell1_pt.Draw("E")
c.SaveAs("h_nu_ell1_pt.pdf")
h_nu_ell2_E.Draw("E")
c.SaveAs("h_nu_ell2_E.pdf")
h_nu_ell2_pt.Draw("E")
c.SaveAs("h_nu_ell2_pt.pdf")

h_mT2_W.Draw("E")
c.SaveAs("h_mT2_W.pdf")

h_mT2_t_11_22.Draw("E")
c.SaveAs("h_mT2_t_11_22.pdf")
h_mT2_t_12_21.Draw("E")
c.SaveAs("h_mT2_t_12_21.pdf")
h_mT2_t_min.Draw("E")
c.SaveAs("h_mT2_t_min.pdf")

h_mT_ell1.Draw("E")
c.SaveAs("h_mT_ell1.pdf")
h_mT_ell1__forMt2Overlay.Draw("E")
c.SaveAs("h_mT_ell1__forMt2Overlay.pdf")

h_EtMiss.Draw("E") 
c.SaveAs("h_EtMiss.pdf")
h_EtMiss_phi.Draw("E") 
c.SaveAs("h_EtMiss_phi.pdf")

h_eventNumber.Draw("E") 
c.SaveAs("h_eventNumber.pdf")

h_truth_numNu.Draw("E")
c.SaveAs("h_truth_numNu.pdf")

h_mT2dc_diff_alpha_1.Draw("E") 
c.SaveAs("h_mT2dc_diff_alpha_1.pdf")
h_mT2dc_alpha_1.Draw("E")
c.SaveAs("h_mT2dc_alpha_1.pdf") 
h_mT2prime_W_alpha_1.Draw("E") 
c.SaveAs("h_mT2prime_W_alpha_1.pdf")
h_mT2prime_t_alpha_1.Draw("E")
c.SaveAs("h_mT2prime_t_alpha_1.pdf")
h_pT_alpha_1.Draw("E")
c.SaveAs("h_pT_alpha_1.pdf")

h_mT2dc_diff_alpha_1_UC.Draw("E") 
c.SaveAs("h_mT2dc_diff_alpha_1_UC.pdf")
h_mT2dc_alpha_1_UC.Draw("E")
c.SaveAs("h_mT2dc_alpha_1_UC.pdf") 
h_mT2prime_W_alpha_1_UC.Draw("E") 
c.SaveAs("h_mT2prime_W_alpha_1_UC.pdf")
h_mT2prime_t_alpha_1_UC.Draw("E")
c.SaveAs("h_mT2prime_t_alpha_1_UC.pdf")
h_pT_sideA_alpha_1_UC.Draw("E")
c.SaveAs("h_pT_sideA_alpha_1_UC.pdf")
h_pT_sideB_alpha_1_UC.Draw("E")
c.SaveAs("h_pT_sideB_alpha_1_UC.pdf")
h_pT_min_alpha_1_UC.Draw("E")
c.SaveAs("h_pT_min_alpha_1_UC.pdf")
h_pT_ratio_alpha_1_UC.Draw("E")
c.SaveAs("h_pT_ratio_alpha_1_UC.pdf")

h_mT2prime_W_alpha_0.Draw("E") 
c.SaveAs("h_mT2prime_W_alpha_0.pdf") 
h_mT2dc_diff_alpha_0.Draw("E") 
c.SaveAs("h_mT2dc_diff_alpha_0.pdf") 
h_mT2dc_alpha_0.Draw("E") 
c.SaveAs("h_mT2dc_alpha_0.pdf") 
h_mT2prime_t_alpha_0.Draw("E") 
c.SaveAs("h_mT2prime_t_alpha_0.pdf")  
h_pT_alpha_0.Draw("E")
c.SaveAs("h_pT_alpha_0.pdf")

h_mT2prime_W_alpha_0_UC.Draw("E") 
c.SaveAs("h_mT2prime_W_alpha_0_UC.pdf") 
h_mT2dc_diff_alpha_0_UC.Draw("E") 
c.SaveAs("h_mT2dc_diff_alpha_0_UC.pdf") 
h_mT2dc_alpha_0_UC.Draw("E") 
c.SaveAs("h_mT2dc_alpha_0_UC.pdf") 
h_mT2prime_t_alpha_0_UC.Draw("E") 
c.SaveAs("h_mT2prime_t_alpha_0_UC.pdf")  
h_pT_sideA_alpha_0_UC.Draw("E")
c.SaveAs("h_pT_sideA_alpha_0_UC.pdf")
h_pT_sideB_alpha_0_UC.Draw("E")
c.SaveAs("h_pT_sideB_alpha_0_UC.pdf")
h_pT_min_alpha_0_UC.Draw("E")
c.SaveAs("h_pT_min_alpha_0_UC.pdf")
h_pT_ratio_alpha_0_UC.Draw("E")
c.SaveAs("h_pT_ratio_alpha_0_UC.pdf")

h_mT2dc_alpha_1_C_vs_UC.Draw("E")
c.SaveAs("h_mT2dc_alpha_1_C_vs_UC.pdf") 

h_mT2dc_alpha_0_C_vs_UC.Draw("E")
c.SaveAs("h_mT2dc_alpha_0_C_vs_UC.pdf") 

h_alpha_1_mT2dcW_vs_minPT.Draw("E")
c.SaveAs("h_alpha_1_mT2dcW_vs_minPT.pdf") 

# save to ROOT output files
h_ell1_pt.Write()
h_ell1_E.Write()

h_ell2_pt.Write()
h_ell2_E.Write()

h_bjet1_E.Write()
h_bjet2_E.Write()

h_nu_ell1_E.Write()
h_nu_ell2_E.Write()

h_nu_ell1_pt.Write()
h_nu_ell2_pt.Write()

h_mT2_W.Write()
h_mT2_t_11_22.Write()
h_mT2_t_12_21.Write()
h_mT2_t_min.Write()

h_mT_ell1.Write()
h_mT_ell1__forMt2Overlay.Write()

h_EtMiss.Write()
h_EtMiss_phi.Write() 

h_eventNumber.Write() 
h_truth_numNu.Write() 

h_mT2dc_diff_alpha_1.Write()
h_mT2dc_alpha_1.Write()
h_mT2prime_W_alpha_1.Write() 
h_mT2prime_t_alpha_1.Write() 
h_pT_alpha_1.Write()

h_mT2dc_diff_alpha_1_UC.Write()
h_mT2dc_alpha_1_UC.Write()
h_mT2prime_W_alpha_1_UC.Write() 
h_mT2prime_t_alpha_1_UC.Write() 
h_pT_sideA_alpha_1_UC.Write()
h_pT_sideB_alpha_1_UC.Write()
h_pT_min_alpha_1_UC.Write() 
h_pT_ratio_alpha_1_UC.Write() 

h_mT2prime_W_alpha_0.Write() 
h_mT2dc_diff_alpha_0.Write()
h_mT2dc_alpha_0.Write() 
h_mT2prime_t_alpha_0.Write() 
h_pT_alpha_0.Write() 

h_mT2prime_W_alpha_0_UC.Write() 
h_mT2dc_diff_alpha_0_UC.Write()
h_mT2dc_alpha_0_UC.Write() 
h_mT2prime_t_alpha_0_UC.Write() 
h_pT_sideA_alpha_0_UC.Write()
h_pT_sideB_alpha_0_UC.Write()
h_pT_min_alpha_0_UC.Write() 
h_pT_ratio_alpha_0_UC.Write() 

h_mT2dc_alpha_1_C_vs_UC.Write() 
h_mT2dc_alpha_1_C_vs_UC.Write() 

h_alpha_1_mT2dcW_vs_minPT.Write()

f_outputRoot.Close()
