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
notable_sol = [] 
notable_sol_2 = [] 

# Define functions
def mass_scalarCalc(px, py, pz, E): 
    """Calculate the mass of a particular, using scalar input variables.
    """
    m_squared = E**2 - (px**2 + py**2 + pz**2)
    
    if m_squared > 0:
        return np.sqrt(m_squared)
    return 0 

def ET_scalarCalc(m, pT):
    """Calculate the transverse energy of a particle, using scalar input variables.
    """
    ET = math.sqrt(max(0, m**2+pT**2)) 
    return ET

def mT_arrayCalc(p4_vis_array, p4_invis_array):
    """Calculate the transverse mass of a parent particle, using array input variables.
    Note: p4_(in)vis_array = [px, py, pz, E].
    """ 
    # Extract mass information
    m_vis = max(0, mass_scalarCalc(p4_vis_array[0], p4_vis_array[1], p4_vis_array[2], p4_vis_array[3]))
    m_invis = max(0, mass_scalarCalc(p4_invis_array[0], p4_invis_array[1], p4_invis_array[2], p4_invis_array[3]))
    
    # Get transverse momentum vectors 
    pT_vis_vec = np.array([p4_vis_array[0], p4_vis_array[1]])  
    pT_invis_vec = np.array([p4_invis_array[0], p4_invis_array[1]])
    
    # Extract energy information 
    ET_vis = ET_scalarCalc(m_invis, np.linalg.norm(pT_vis_vec)) 
    ET_invis = ET_scalarCalc(m_invis, np.linalg.norm(pT_invis_vec)) 
   
    mT = math.sqrt(max(0, m_vis**2 + m_invis**2 + 2*(ET_vis*ET_invis - np.dot(pT_vis_vec, pT_invis_vec))))
    return mT 

# Define TLorentzVector module functions 
def extract_Px(pT, eta, phi, mass): 
    """Extract Px, from scalar input variables. 
    """
    Px = pT*np.cos(phi) 
    return Px 

def extract_Py(pT, eta, phi, mass):
    """Extract Py, from scalar input variables.
    """
    Py = pT*np.sin(phi) 
    return Py 
    
def extract_Pz(pT, eta, phi, mass):
    """Extract Pz, from scalar input variables.
    """ 
    theta = 2*np.arctan(np.exp(-eta)) # polar angle 
    Pz = pT/np.tan(theta)
    return Pz 
    
def extract_E(pT, eta, phi, mass):
    """Extract E, from scalar input variables.
    """
    Pz = extract_Pz(pT, eta, phi, mass)
    E = math.sqrt(mass**2 + (pT**2 + Pz**2))
    return E 
    
# Create TH1 histogram 
# mT2dc(alpha = 1) - mT2(W)   
h_alpha_1 = ROOT.TH1F("h_alpha_1", "mT2dc(alpha = 1) - mt2dc; Difference [GeV]; Number of entries / 1 GeV", 100, -100, 100)
h_alpha_2 = ROOT.TH1F("h_alpha_2", "mT2dc(alpha = 1) - mt2dc; Difference [GeV]; Number of entries / 1 GeV", 100, -100, 100)
h_alpha_3 = ROOT.TH1F("h_alpha_3", "mT2dc(alpha = 1) - mt2dc; Difference [GeV]; Number of entries / 1 GeV", 100, -100, 100)
h_alpha_4 = ROOT.TH1F("h_alpha_4", "mT2dc(alpha = 1) - mt2dc; Difference [GeV]; Number of entries / 1 GeV", 100, -100, 100)
h_alpha_5 = ROOT.TH1F("h_alpha_5", "mT2dc(alpha = 1) - mt2dc; Difference [GeV]; Number of entries / 1 GeV", 100, -100, 100)
h_alpha_6 = ROOT.TH1F("h_alpha_6", "mT2dc(alpha = 1) - mt2dc; Difference [GeV]; Number of entries / 1 GeV", 100, -100, 100)
h_alpha_7 = ROOT.TH1F("h_alpha_7", "mT2dc(alpha = 1) - mt2dc; Difference [GeV]; Number of entries / 1 GeV", 100, -100, 100)

h_mt2dc_sol_1 = ROOT.TH1F("h_m2tdc_sol_1", "mT2dc(alpha = 1); mT2dc [GeV]; Number of entries / 1 GeV", 100, 0, 100)
h_mt2dc_sol_2 = ROOT.TH1F("h_m2tdc_sol_2", "mT2dc(alpha = 1); mT2dc [GeV]; Number of entries / 1 GeV", 100, 0, 100)
h_mt2dc_sol_3 = ROOT.TH1F("h_m2tdc_sol_3", "mT2dc(alpha = 1); mT2dc [GeV]; Number of entries / 1 GeV", 100, 0, 100)
h_mt2dc_sol_4 = ROOT.TH1F("h_m2tdc_sol_4", "mT2dc(alpha = 1); mT2dc [GeV]; Number of entries / 1 GeV", 100, 0, 100)
h_mt2dc_sol_5 = ROOT.TH1F("h_m2tdc_sol_5", "mT2dc(alpha = 1); mT2dc [GeV]; Number of entries / 1 GeV", 100, 0, 100)
h_mt2dc_sol_6 = ROOT.TH1F("h_m2tdc_sol_6", "mT2dc(alpha = 1); mT2dc [GeV]; Number of entries / 1 GeV", 100, 0, 100)
h_mt2dc_sol_7 = ROOT.TH1F("h_m2tdc_sol_7", "mT2dc(alpha = 1); mT2dc [GeV]; Number of entries / 1 GeV", 100, 0, 100)

# Get the number entries in the tree 
nentries = t.GetEntries() # 60599  

for i in range(10000):
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
    bjet1_sideA_Px = extract_Px(t.bjet1_PT, t.bjet1_Eta, t.bjet1_Phi, t.bjet1_Mass) 
    bjet1_sideA_Py = extract_Py(t.bjet1_PT, t.bjet1_Eta, t.bjet1_Phi, t.bjet1_Mass)
    bjet1_sideA_Pz = extract_Pz(t.bjet1_PT, t.bjet1_Eta, t.bjet1_Phi, t.bjet1_Mass)
    bjet1_sideA_E = extract_E(t.bjet1_PT, t.bjet1_Eta, t.bjet1_Phi, t.bjet1_Mass) 
    bjet1_sideA_array = np.array([bjet1_sideA_Px, bjet1_sideA_Py, bjet1_sideA_Pz, bjet1_sideA_E]) 
    
    # get sideA lepton information 
    ell1_sideA_Px = extract_Px(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0) 
    ell1_sideA_Py = extract_Py(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0)
    ell1_sideA_Pz = extract_Pz(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0)
    ell1_sideA_E = extract_E(t.ell1_PT, t.ell1_Eta, t.ell1_Phi, 0) 
    ell1_sideA_array = np.array([ell1_sideA_Px, ell1_sideA_Py, ell1_sideA_Pz, ell1_sideA_E]) 
    
    vis_sideA_array = np.array([bjet1_sideA_array, ell1_sideA_array]) 
    
    # get sideB bjet information
    bjet2_sideB_Px = extract_Px(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass) 
    bjet2_sideB_Py = extract_Py(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass)
    bjet2_sideB_Pz = extract_Pz(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass)
    bjet2_sideB_E = extract_E(t.bjet2_PT, t.bjet2_Eta, t.bjet2_Phi, t.bjet2_Mass) 
    bjet2_sideB_array = np.array([bjet2_sideB_Px, bjet2_sideB_Py, bjet2_sideB_Pz, bjet2_sideB_E]) 
    
    # get sideB lepton information 
    ell2_sideB_Px = extract_Px(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0) 
    ell2_sideB_Py = extract_Py(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0)
    ell2_sideB_Pz = extract_Pz(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0)
    ell2_sideB_E = extract_E(t.ell2_PT, t.ell2_Eta, t.ell2_Phi, 0) 
    ell2_sideB_array = np.array([ell2_sideB_Px, ell2_sideB_Py, ell2_sideB_Pz, ell2_sideB_E])  
   
    vis_sideB_array = np.array([bjet2_sideB_array, ell2_sideB_array]) 

    # get met information 
    met_Px = extract_Px(t.EtMiss, 0, t.EtMiss_phi, 0) 
    met_Py = extract_Py(t.EtMiss, 0, t.EtMiss_phi, 0) 
    met_E = extract_E(t.EtMiss, 0, t.EtMiss_phi, 0) 
    met = np.array([met_Px, met_Py, 0, met_E])
   
    # define initial solution vector 
    invis_sideA_array_guess = met/2
    invis_sideA_array_guess_2 = bjet1_sideA_array
    invis_sideA_array_guess_3 = ell1_sideA_array 
    invis_sideA_array_guess_4 = bjet2_sideB_array
    invis_sideA_array_guess_5 = ell2_sideB_array

    truth_nu_ell1_Px = extract_Px(t.truth_nu_ell1_PT, t.truth_nu_ell1_Eta, t.truth_nu_ell1_Phi, 0) 
    truth_nu_ell1_Py = extract_Py(t.truth_nu_ell1_PT, t.truth_nu_ell1_Eta, t.truth_nu_ell1_Phi, 0)  
    truth_nu_ell1_Pz = extract_Pz(t.truth_nu_ell1_PT, t.truth_nu_ell1_Eta, t.truth_nu_ell1_Phi, 0) 
    truth_nu_ell1_E = extract_E(t.truth_nu_ell1_PT, t.truth_nu_ell1_Eta, t.truth_nu_ell1_Phi, 0) 
    invis_sideA_array_guess_6 = np.array([truth_nu_ell1_Px, truth_nu_ell1_Py, truth_nu_ell1_Pz, truth_nu_ell1_E]) 
    
    truth_nu_ell2_Px = extract_Px(t.truth_nu_ell2_PT, t.truth_nu_ell2_Eta, t.truth_nu_ell2_Phi, 0) 
    truth_nu_ell2_Py = extract_Py(t.truth_nu_ell2_PT, t.truth_nu_ell2_Eta, t.truth_nu_ell2_Phi, 0)  
    truth_nu_ell2_Pz = extract_Pz(t.truth_nu_ell1_PT, t.truth_nu_ell2_Eta, t.truth_nu_ell2_Phi, 0) 
    truth_nu_ell2_E = extract_E(t.truth_nu_ell2_PT, t.truth_nu_ell2_Eta, t.truth_nu_ell2_Phi, 0) 
    invis_sideA_array_guess_7 = np.array([truth_nu_ell2_Px, truth_nu_ell2_Py, truth_nu_ell2_Pz, truth_nu_ell2_E]) 
    
    # define the function to minimise; minimise over two variables 
    def objective(invis_sideA_array):
        alpha_term_1 = mT_arrayCalc(vis_sideA_array[-1], invis_sideA_array) # mT(lA, pT_A)
        alpha_term_2 = mT_arrayCalc(vis_sideB_array[-1], met-invis_sideA_array) # mT(TB, pT_B) 
        alpha_term = max(alpha_term_1, alpha_term_2) 
        
        beta_term_1 = mT_arrayCalc(vis_sideA_array[0] + vis_sideA_array[-1], invis_sideA_array) # mT(lATA, pT_A)
        beta_term_2 = mT_arrayCalc(vis_sideB_array[0] + vis_sideB_array[-1], met-invis_sideA_array) # mT(TBbB, pt_B)
        beta_term = max(beta_term_1, beta_term_2) 
    
        return alphaList[0]*alpha_term + (1-alphaList[0])*beta_term 

    sol_1 = so.minimize(objective, x0 = invis_sideA_array_guess, method='Nelder-Mead', 
                            options={'maxiter': 2000, 'xatol': 1e-5, 'fatol': 1e-5, 'adaptive': True, 'disp': True}) 
   
    sol_2 = so.minimize(objective, x0 = invis_sideA_array_guess_2, method='Nelder-Mead', 
                            options={'maxiter': 2000, 'xatol': 1e-5, 'fatol': 1e-5, 'adaptive': True})  
    
    sol_3 = so.minimize(objective, x0 = invis_sideA_array_guess_3, method='Nelder-Mead', 
                            options={'maxiter': 2000, 'xatol': 1e-5, 'fatol': 1e-5, 'adaptive': True})      
    
    sol_4 = so.minimize(objective, x0 = invis_sideA_array_guess_4, method='Nelder-Mead', 
                            options={'maxiter': 2000, 'xatol': 1e-5, 'fatol': 1e-5, 'adaptive': True})  
  
    sol_5 = so.minimize(objective, x0 = invis_sideA_array_guess_5, method='Nelder-Mead', 
                            options={'maxiter': 2000, 'xatol': 1e-5, 'fatol': 1e-5, 'adaptive': True})  
    
    sol_6 = so.minimize(objective, x0 = invis_sideA_array_guess_6, method='Nelder-Mead', 
                            options={'maxiter': 2000, 'xatol': 1e-5, 'fatol': 1e-5, 'adaptive': True}) 

    sol_7 = so.minimize(objective, x0 = invis_sideA_array_guess_7, method='Nelder-Mead', 
                            options={'maxiter': 2000, 'xatol': 1e-5, 'fatol': 1e-5, 'adaptive': True})  
   
    print('1', 'event', i, sol_1.fun - mt2_W) 
    print('2', 'event', i, sol_2.fun - mt2_W) 
    print('3', 'event', i, sol_3.fun - mt2_W) 
    print('4', 'event', i, sol_4.fun - mt2_W) 
    print('5', 'event', i, sol_5.fun - mt2_W)     
    print('6', 'event', i, sol_6.fun - mt2_W) 
    print('7', 'event', i, sol_7.fun - mt2_W) 
    
    h_alpha_1.Fill(sol_1.fun - mt2_W) 
    h_mt2dc_sol_1.Fill(sol_1.fun)

    h_alpha_2.Fill(sol_2.fun - mt2_W) 
    h_mt2dc_sol_2.Fill(sol_2.fun)

    h_alpha_3.Fill(sol_3.fun - mt2_W) 
    h_mt2dc_sol_3.Fill(sol_3.fun)

    h_alpha_4.Fill(sol_4.fun - mt2_W) 
    h_mt2dc_sol_4.Fill(sol_4.fun)
    
    h_alpha_5.Fill(sol_5.fun - mt2_W) 
    h_mt2dc_sol_5.Fill(sol_5.fun)
    
    h_alpha_6.Fill(sol_6.fun - mt2_W) 
    h_mt2dc_sol_6.Fill(sol_6.fun)
       
    h_alpha_7.Fill(sol_7.fun - mt2_W) 
    h_mt2dc_sol_7.Fill(sol_7.fun) 
    
    #if sol.fun - mt2_W > 1:
        #notable_sol.append([i, sol.fun - mt2_W])
    #if sol_2.fun - mt2_W > 1:
       # notable_sol_2.append([i, sol_2.fun - mt2_W]) 
        
    # Get alpha, beta terms 
    #def get_alpha_term():
    #    alpha_term_1 = mT_arrayCalc(vis_sideA_array[-1], sol.x) # mT(lA, pT_A)
    #    alpha_term_2 = mT_arrayCalc(vis_sideB_array[-1], met-sol.x) # mT(TB, pT_B) 
    #    alpha_term = max(alpha_term_1, alpha_term_2)
    #    return alpha_term 
   
    #def get_beta_term():
    #    beta_term_1 = mT_arrayCalc(vis_sideA_array[0] + vis_sideA_array[-1], sol.x) # mT(lATA, pT_A)
    #    beta_term_2 = mT_arrayCalc(vis_sideB_array[0] + vis_sideB_array[-1], met-sol.x) # mT(TBbB, pt_B)
     #   beta_term = max(beta_term_1, beta_term_2) 
    #    return beta_term 
    
    #def get_invisible_directions(): 
    #    invis_sideA_2vec = np.array([sol.x[0], sol.x[1]])  
    #    invis_sideB_2vec = np.array([met[0], met[1]]) - np.array([sol.x[0], sol.x[1]]) 
    #    return invis_sideA_2vec, invis_sideB_2vec 
    
    #print('alpha_term', get_alpha_term()) 
    #print('beta_term', get_beta_term()) 
    #print('invis_sideA, invis_sideB directions', get_invisible_directions()) 
                                  
#print(notable_sol) 

# Draw the histograms and save them.
c = ROOT.TCanvas()
                         
h_alpha_1.Draw("E") # put error bars 
c.SaveAs("h_alpha_1_met_METHOD1_e2.pdf")

h_mt2dc_sol_1.Draw("E") # put error bars 
c.SaveAs("h_mt2dc_sol_met_METHOD1.pdf")

h_alpha_2.Draw("E") # put error bars 
c.SaveAs("h_alpha_1_met_METHOD2.pdf")

h_mt2dc_sol_2.Draw("E") # put error bars 
c.SaveAs("h_mt2dc_sol_met_METHOD2.pdf")

h_alpha_3.Draw("E") # put error bars 
c.SaveAs("h_alpha_1_met_METHOD3.pdf")

h_mt2dc_sol_3.Draw("E") # put error bars 
c.SaveAs("h_mt2dc_sol_met_METHOD3.pdf")

h_alpha_4.Draw("E") # put error bars 
c.SaveAs("h_alpha_1_met_METHOD4.pdf")

h_mt2dc_sol_4.Draw("E") # put error bars 
c.SaveAs("h_mt2dc_sol_met_METHOD4.pdf")

h_alpha_5.Draw("E") # put error bars 
c.SaveAs("h_alpha_1_met_METHOD5.pdf")

h_mt2dc_sol_5.Draw("E") # put error bars 
c.SaveAs("h_mt2dc_sol_met_METHOD5.pdf")

h_alpha_6.Draw("E") # put error bars 
c.SaveAs("h_alpha_1_met_METHOD6.pdf")

h_mt2dc_sol_6.Draw("E") # put error bars 
c.SaveAs("h_mt2dc_sol_met_METHOD6.pdf")


h_alpha_7.Draw("E") # put error bars 
c.SaveAs("h_alpha_1_met_METHOD7.pdf")

h_mt2dc_sol_7.Draw("E") # put error bars 
c.SaveAs("h_mt2dc_sol_met_METHOD7.pdf")


h_alpha_1.Write() 
h_mt2dc_sol_1.Write() 

f_outputRoot.Close()
