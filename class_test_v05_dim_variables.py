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

h_mt2dc_sol = ROOT.TH1F("h_m2tdc_sol", "mT2dc(alpha = 1); mT2dc [GeV]; Number of entries / 1 GeV", 200, 0, 200)

# Get the number entries in the tree 
nentries = t.GetEntries() # 60599  

for i in range(5000):
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
   
    # define initial solution vector guesses 
    invis_sideA_array_guess = [met[0], met[1], 0, 0]/2 
    invis_sideA_array_guess_2 = [bjet1_sideA_array[0], bjet1_sideA_array[1], 0, 0] 
    invis_sideA_array_guess_3 = ell1_sideA_array[:2] 
    invis_sideA_array_guess_4 = bjet2_sideB_array[:2] 
    invis_sideA_array_guess_5 = ell2_sideB_array[:2] 
    invis_sideA_array_guess_6 = 1.01*met[:2] 
    invis_sideA_array_guess_7 = 0.01*met[:2] 
    
    invis_sideA_array_guesses = [invis_sideA_array_guess, invis_sideA_array_guess_2, invis_sideA_array_guess_3, 
                                invis_sideA_array_guess_4, invis_sideA_array_guess_5, invis_sideA_array_guess_6, 
                                invis_sideA_array_guess_7]

    # define the function to minimise; minimise over two variables 
    def objective(invis_sideA_2vec, dim_variable_1, dim_variable_2): # minimise over a 2-vector array, having components px and py 
        invis_sideA_array = np.array([invis_sideA_2vec[0], invis_sideA_2vec[1], 0, 
                                     np.sqrt(invis_sideA_2vec[0]**2 + invis_sideA_2vec[0]**2)]) 
        
        alpha_term_1 = mT_arrayCalc(vis_sideA_array[-1], invis_sideA_array) # mT(lA, pT_A)
        alpha_term_2 = mT_arrayCalc(vis_sideB_array[-1], met-invis_sideA_array) # mT(TB, pT_B) 
        alpha_term = max(alpha_term_1, alpha_term_2) 
        
        beta_term_1 = mT_arrayCalc(vis_sideA_array[0] + vis_sideA_array[-1], invis_sideA_array) # mT(lATA, pT_A)
        beta_term_2 = mT_arrayCalc(vis_sideB_array[0] + vis_sideB_array[-1], met-invis_sideA_array) # mT(TBbB, pt_B)
        beta_term = max(beta_term_1, beta_term_2) 
    
        return alphaList[0]*alpha_term + (1-alphaList[0])*beta_term 
    
    # try different inputs to establish which produces the smallest value
    guess = objective(invis_sideA_array_guess)
    guess_2 = objective(invis_sideA_array_guess_2) 
    guess_3 = objective(invis_sideA_array_guess_3) 
    guess_4 = objective(invis_sideA_array_guess_4)
    guess_5 = objective(invis_sideA_array_guess_5) 
    guess_6 = objective(invis_sideA_array_guess_6)
    guess_7 = objective(invis_sideA_array_guess_7)
    
    guesses = [guess, guess_2, guess_3, guess_4, guess_5, guess_6, guess_7] 
    #print('guesses', [guess_1, guess_2, guess_3, guess_4, guess_5, guess_6, guess_7]) 
    #print('argmin', np.argmin(guesses), guesses[np.argmin(guesses)], invis_sideA_array_guesses[np.argmin(guesses)]) 
    
    sol = so.minimize(objective, x0 = invis_sideA_array_guesses[np.argmin(guesses)], method='Nelder-Mead', 
                            options={'maxiter': 2500, 'xatol': 2e-6, 'fatol': 2e-6, 'adaptive': True, 'disp': True}) 
   
    print('event', i, sol.fun - mt2_W) 
    
    h_alpha_1.Fill(sol.fun - mt2_W) 
    h_mt2dc_sol.Fill(sol.fun)
    
#Draw the histograms and save them.
c = ROOT.TCanvas()
                         
h_alpha_1.Draw("E") # put error bars 
c.SaveAs("h_alpha_1_NM_2vec_increase.pdf")

h_mt2dc_sol.Draw("E") # put error bars 
c.SaveAs("h_mt2dc_sol_NM_2vec_increase.pdf")

h_alpha_1.Write() 
h_mt2dc_sol.Write() 

f_outputRoot.Close()
