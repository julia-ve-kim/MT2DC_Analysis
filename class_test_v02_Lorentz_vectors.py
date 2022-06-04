#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 10:02:03 2022

@author: juliakim
"""

import ROOT
import math
import numpy as np
import scipy.optimize as so 

# Define input file 
f_inputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_03_March_07_skim_mg5_ttbar_jet_merged_001-716_ntuple_2l2b_v01.root", "read")
t = f_inputRoot.Get("variables")
type(t)
f_outputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_05_May_10_mt2dc_analysis_v01.root", "recreate")

# Define constants 
alphaList = [1] 

# Define functions
def Et_scalarCalc(m, pt):
    """Calculate the transverse energy of a particle, using the scalar input variables: mass (m) and transverse momentum (pt)."""
    Et = math.sqrt(max(0, m**2+pt**2)) 
    return Et

# Operations Required for TLorentz Vectors 
def mT_4vecCalc(p4_vis_array, p4_invis_array):
    """Calculate the tranverse mass of a parent particle, using array input variables. 
    """
  
    # Create TLorentz 4-vectors 
    p4_vis = ROOT.TLorentzVector() 
    p4_vis.SetPxPyPzE(p4_vis_array[0], p4_vis_array[1], p4_vis_array[2], p4_vis_array[3])
   
    p4_invis = ROOT.TLorentzVector() 
    p4_invis.SetPxPyPzE(p4_invis_array[0], p4_invis_array[1], p4_invis_array[2], p4_invis_array[3])
    
    # Extract mass, energy and transverse momentum information 
    m_vis = max(0, p4_vis.M()) 
    m_invis = max(0, p4_invis.M()) 
    
    Et_vis = Et_scalarCalc(m_vis, p4_vis.Pt()) 
    Et_invis = Et_scalarCalc(m_invis, p4_invis.Pt()) 
    
    pt_vec_vis = p4_vis.Vect() 
    pt_vec_vis.SetZ(0)
    pt_vec_invis = p4_invis.Vect() 
    pt_vec_invis.SetZ(0) 
    
    mT = math.sqrt(max(0, m_vis**2 + m_invis**2 + 2*(Et_vis*Et_invis - pt_vec_vis*pt_vec_invis)))
    return mT 

# Create a TH1 histogram 
# mT2(W) - mT2dc(alpha = 0)  
h_alpha_1_v02 = ROOT.TH1F("h_alpha_1", "mT2dc(alpha = 1) - mt2DC; Difference (GeV); Number of entries / 1 GeV", 100, -100, 100)
h_sol_array = ROOT.TH1F("h_sol_array", "mT2dc(alpha = 1); mT2dc (GeV); Number of entries / 1 GeV", 100, 0, 100)
notable_sol = [] 
    
# Get the number entries in the tree 
nentries = t.GetEntries() # 60599  

# Testing Input #2 
for i in range(0, 1000):
    if (i%1000==0): 
       print(":: Processing entry ", i, " = ", i*1.0/nentries*100.0, "%.")    
    if t.LoadTree(i) < 0:
       print("**could not load tree for entry #%s") % i
       break
    nb = t.GetEntry(i) 
    if nb <= 0:
       continue
    
    # retrive information from tree 
    # get the mt2_W information from every event 
    mt2_W = t.mt2_W_ell1ell2 
    
    # get sideA bjet information 
    bjet1_sideA_p4 = ROOT.TLorentzVector() 
    bjet1_sideA_p4.SetPtEtaPhiM(t.bjet1_PT, t.bjet1_Eta, t.bjet1_Phi, t.bjet1_Mass) 
    bjet1_sideA_array = np.array([bjet1_sideA_p4.Px(), bjet1_sideA_p4.Py(), bjet1_sideA_p4.Pz(), bjet1_sideA_p4.E()]) 
    
    # get sideA lepton information 
    ell1_sideA_p4 = ROOT.TLorentzVector() 
    ell1_sideA_p4.SetPtEtaPhiM(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0)
    ell1_sideA_array = np.array([ell1_sideA_p4.Px(), ell1_sideA_p4.Py(), ell1_sideA_p4.Pz(), ell1_sideA_p4.E()])
    
    vis_sideA_array = np.array([bjet1_sideA_array, ell1_sideA_array]) 
    
    # get sideB bjet information
    bjet2_sideB_p4 = ROOT.TLorentzVector() 
    bjet2_sideB_p4.SetPtEtaPhiM(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass) 
    bjet2_sideB_array = np.array([bjet2_sideB_p4.Px(), bjet2_sideB_p4.Py(), bjet2_sideB_p4.Pz(), bjet2_sideB_p4.E()])
    
    # get sideB lepton information 
    ell2_sideB_p4 = ROOT.TLorentzVector() 
    ell2_sideB_p4.SetPtEtaPhiM(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0) 
    ell2_sideB_array = np.array([ell2_sideB_p4.Px(), ell2_sideB_p4.Py(), ell2_sideB_p4.Pz(), ell2_sideB_p4.E()])
    
    vis_sideB_array = np.array([bjet2_sideB_array, ell2_sideB_array]) 

    # get met information 
    met_p4 = ROOT.TLorentzVector()
    met_p4.SetPtEtaPhiM(t.EtMiss, 0, t.EtMiss_phi, 0)
    met = np.array([met_p4.Px(), met_p4.Py(), 0, met_p4.E()]) 
    
    # define initial solution vector 
    invis_sideA_array_guess = met/2 
    
    # define the function to minimise
    def objective(invis_sideA_array):
        alpha_term_1 = mT_4vecCalc(vis_sideA_array[-1], invis_sideA_array) # mT(lA, pT_A)
        alpha_term_2 = mT_4vecCalc(vis_sideB_array[-1], met-invis_sideA_array) # mT(TB, pT_B) 
        alpha_term = max(alpha_term_1, alpha_term_2) 
        
        beta_term_1 = mT_4vecCalc(vis_sideA_array[0] + vis_sideA_array[-1], invis_sideA_array) # mT(lATA, pT_A)
        beta_term_2 = mT_4vecCalc(vis_sideB_array[0] + vis_sideB_array[-1], met-invis_sideA_array) # mT(TBbB, pt_B)
        beta_term = max(beta_term_1, beta_term_2) 
    
        return alphaList[0]*alpha_term + (1-alphaList[0])*beta_term 

    sol_array_1 = so.minimize(objective, x0 = invis_sideA_array_guess, method='Nelder-Mead', 
                            options={'xatol': 1e-15, 'fatol': 1e-15, 'adaptive': True})
    
    sol_array_2 = so.minimize(objective, x0 = invis_sideA_array_guess, method='Nelder-Mead', 
                            options={'xatol': 1e-15, 'fatol': 1e-15, 'adaptive': False})
    
    sol_fun = min(sol_array_1.fun, sol_array_2.fun) 
    
    print('event', i, sol_fun - mt2_W) # event numbers + print mt2_w - sol_array.fun
    # fill histogram 
    h_alpha_1_v02.Fill(sol_fun - mt2_W) 
    h_sol_array.Fill(sol_fun)
    
    if sol_fun - mt2_W > 1:
        notable_sol.append([i, sol_fun - mt2_W])
                         
# Draw the histograms and save them.
print(notable_sol) 
c = ROOT.TCanvas()
                         
h_alpha_1_v02.Draw("E") # put error bars 
c.SaveAs("h_alpha_1_met_refined_test_3.pdf")

h_sol_array.Draw("E") # put error bars 
c.SaveAs("h_alpha_1_sol_refined_test_3.pdf")

h_alpha_1_v02.Write() 
h_sol_array.Write() 

f_outputRoot.Close()
