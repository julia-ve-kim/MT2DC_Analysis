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
f_outputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_05_May_10_mt2dc_analysis_v01.root", "RECREATE")

##############################################
# Define the plots to produce
##############################################
# Plots from input ROOT TFile Tree                 

h_mT2_W = ROOT.TH1F("h_mT2_W", "mt2(ell1,ell2) = mt2(W); mt2(W) [GeV]; Number of entries / 1 GeV", 200, 0, 200)
h_mT2_t_11_22 = ROOT.TH1F("h_mT2_t_11_22", "mt2|11,22(t) = mt2(b1 ell1,b2 ell2); mt2(t|11,22) [GeV]; Number of entries / 1 GeV", 300, 0, 300)

h_mT2dc_diff_alpha_1 = ROOT.TH1F("h_mT2dc_diff_alpha_1", "mT2dc(alpha = 1) - mT2(W); Difference [GeV]; Number of entries / 2 GeV", 100, -100, 100)
h_mT2dc_alpha_1 = ROOT.TH1F("h_mT2dc_alpha_1", "mT2dc(alpha = 1); mT2dc [GeV]; Number of entries / 1 GeV", 200, 0, 200)
h_mT2prime_W_alpha_1 = ROOT.TH1F("h_mT2prime_W_alpha_1", "mt2'(W); mt2'(W) [GeV]; Number of entries / 1 GeV", 200, 0, 200)
h_mT2prime_t_alpha_1 = ROOT.TH1F("h_mT2prime_t_alpha_1", "mt2'(W); mt2'(t) [GeV]; Number of entries / 1 GeV", 200, 0, 200)

h_mT2dc_diff_alpha_0  = ROOT.TH1F("h_mT2dc_diff_alpha_0", "mT2dc(alpha = 0) - mt2_t_bjet1ell1_bjet2ell2; Difference [GeV]; Number of entries / 2 GeV", 100, -100, 100)
h_mT2dc_alpha_0 = ROOT.TH1F("h_mT2dc_alpha_0", "mT2dc(alpha = 0); mT2dc [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_mT2prime_W_alpha_0 = ROOT.TH1F("h_mT2prime_W_alpha_0", "mt2'(W); mt2'(W) [GeV]; Number of entries / 1 GeV", 300, 0, 300)
h_mT2prime_t_alpha_0 = ROOT.TH1F("h_mT2prime_t", "mt2'(t); mt2'(t) [GeV]; Number of entries / 1 GeV", 300, 0, 300)

##############################################
# Define constants
##############################################
m_W = 80.   # GeV 
m_t = 173.  # GeV
nentries = t.GetEntries() # 60599
alphaListindex = input('Choose 0 (alpha = 0) or 1 (alpha = 1):') 
calcStyle = input('Choose fast or slow:') 

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
    
    # get sideA lepton information (suppose sideA = highest p light lepton) 
    ell1_sideA_Px = DC.extract_Px(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0) 
    ell1_sideA_Py = DC.extract_Py(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0)
    ell1_sideA_Pz = DC.extract_Pz(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0)
    ell1_sideA_E = DC.extract_E(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0) 
    ell1_sideA_array = np.array([ell1_sideA_Px, ell1_sideA_Py, ell1_sideA_Pz, ell1_sideA_E])
    
    vis_sideA_array = np.array([bjet1_sideA_array, ell1_sideA_array]) 
    
    # get sideB bjet information
    bjet2_sideB_Px = DC.extract_Px(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass) 
    bjet2_sideB_Py = DC.extract_Py(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass)
    bjet2_sideB_Pz = DC.extract_Pz(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass)
    bjet2_sideB_E = DC.extract_E(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass) 
    bjet2_sideB_array = np.array([bjet2_sideB_Px, bjet2_sideB_Py, bjet2_sideB_Pz, bjet2_sideB_E])
    
    # get sideB lepton information (suppose sideB = second highest p light lepton) 
    ell2_sideB_Px = DC.extract_Px(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0) 
    ell2_sideB_Py = DC.extract_Py(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0)
    ell2_sideB_Pz = DC.extract_Pz(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0)
    ell2_sideB_E = DC.extract_E(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0) 
    ell2_sideB_array = np.array([ell2_sideB_Px, ell2_sideB_Py, ell2_sideB_Pz, ell2_sideB_E]) 
        
    vis_sideB_array = np.array([bjet2_sideB_array, ell2_sideB_array]) 
    
    # get met information 
    met_Px = DC.extract_Px(t.EtMiss, 0, t.EtMiss_phi, 0) 
    met_Py = DC.extract_Py(t.EtMiss, 0, t.EtMiss_phi, 0) 
    met_E = DC.extract_E(t.EtMiss, 0, t.EtMiss_phi, 0) 
    met = np.array([met_Px, met_Py, 0, met_E])

    # get tranverse mass information
    mt2_W = t.mt2_W_ell1ell2
    mt2_t_11_22 = t.mt2_t_bjet1ell1_bjet2ell2

    h_mT2_W.Fill(mt2_W)
    h_mT2_t_11_22.Fill(mt2_t_11_22)
   
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
    
    if calcStyle == 'fast': 
        invis_sideA_array_guesses = [invis_sideA_array_guess, invis_sideA_array_guess_2, invis_sideA_array_guess_3, 
                                     invis_sideA_array_guess_4, invis_sideA_array_guess_5, invis_sideA_array_guess_6, 
                                     invis_sideA_array_guess_7]
        guess = objective(invis_sideA_array_guess)
        guess_2 = objective(invis_sideA_array_guess_2) 
        guess_3 = objective(invis_sideA_array_guess_3) 
        guess_4 = objective(invis_sideA_array_guess_4)
        guess_5 = objective(invis_sideA_array_guess_5) 
        guess_6 = objective(invis_sideA_array_guess_6)
        guess_7 = objective(invis_sideA_array_guess_7)
        
        guesses = [guess, guess_2, guess_3, guess_4, guess_5, guess_6, guess_7] 
        sol = so.minimize(objective, x0 = invis_sideA_array_guesses[np.argmin(guesses)], method='SLSQP', 
                          options={'maxiter': 2000, 'ftol': 1e-07,'disp': True}) 
        
        if index == 0:
            h_mT2dc_diff_alpha_0.Fill(sol_fun - mt2_t_11_22) 
            h_mT2dc_alpha_0.Fill(sol_fun) 
            h_mT2prime_W_alpha_0.Fill(DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_x)) 
            h_mT2prime_t_alpha_0.Fill(DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_x))  
        elif index == 1:
            h_mT2dc_diff_alpha_1.Fill(sol_fun - mt2_W) 
            h_mT2dc_alpha_1.Fill(sol_fun) 
            h_mT2prime_W_alpha_1.Fill(DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_x))  
            h_mT2prime_t_alpha_1.Fill(DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_x))  
    
    elif calcStyle == 'slow':
        sol_1 = so.minimize(objective, x0 = invis_sideA_array_guess, method='SLSQP', 
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
            h_mT2prime_W_alpha_0.Fill(DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_x)) 
            h_mT2prime_t_alpha_0.Fill(DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_x))  
        elif index == 1:
            h_mT2dc_diff_alpha_1.Fill(sol_fun - mt2_W) 
            h_mT2dc_alpha_1.Fill(sol_fun) 
            h_mT2prime_W_alpha_1.Fill(DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_x))  
            h_mT2prime_t_alpha_1.Fill(DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_x))  
            
            
##############################################
# Draw all histograms and save them.
##############################################

c = ROOT.TCanvas()

h_mT2_W.Draw("E")
c.SaveAs("h_mT2_W.pdf")
h_mT2_t_11_22.Draw("E")
c.SaveAs("h_mT2_t_11_22.pdf")

h_mT2dc_diff_alpha_1.Draw("E") 
c.SaveAs("h_mT2dc_diff_alpha_1.pdf")
h_mT2dc_alpha_1.Draw("E")
c.SaveAs("h_mT2dc_alpha_1.pdf") 
h_mT2prime_W_alpha_1.Draw("E")
c.SaveAs("h_mT2prime_W_alpha_1.pdf")
h_mT2prime_t_alpha_1.Draw("E")
c.SaveAs("h_mT2prime_t_alpha_1.pdf")

h_mT2dc_diff_alpha_0.Draw("E") 
c.SaveAs("h_mT2dc_diff_alpha_0.pdf") 
h_mT2dc_alpha_0.Draw("E") 
c.SaveAs("h_mT2dc_alpha_0.pdf") 
h_mT2prime_t_alpha_0.Draw("E") 
c.SaveAs("h_mT2prime_t_alpha_0.pdf") 
h_mT2prime_W_alpha_0.Draw("E") 
c.SaveAs("h_mT2prime_W_alpha_0.pdf") 
    
# save to ROOT output files
h_mT2_W.Write()
h_mT2_t_11_22.Write()

h_mT2dc_diff_alpha_1.Write()
h_mT2dc_alpha_1.Write()
h_mT2prime_W_alpha_1.Write() 
h_mT2prime_t_alpha_1.Write() 

h_mT2dc_diff_alpha_0.Write()
h_mT2dc_alpha_0.Write() 
h_mT2prime_W_alpha_0.Write() 
h_mT2prime_t_alpha_0.Write() 

f_outputRoot.Close()
