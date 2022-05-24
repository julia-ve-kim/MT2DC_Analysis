#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 10:02:03 2022

@author: juliakim
"""

import ROOT
import mt2dc as DC
import math
import numpy as np
import scipy.optimize as so 

# Class Attributes  
def Et_scalarCalc(m, pt):
    """Calculate the transverse energy of a particle, using scalar input variables."""
    Et = math.sqrt(max(0, m**2+pt**2)) 
    return Et

def mT_4vecCalc(p4_vis_array, p4_invis_array):
    """Calculate the tranverse mass of a parent particle, using array input variables. 
    Inputs: p4_vis_array = [px, py, pz, E], p4_invis_array = [px, py, 0, E]. 
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
   #pt_vec_invis.SetZ(0) 
    
    mT = math.sqrt(max(0, m_vis**2 + m_invis**2 + 2*(Et_vis*Et_invis - pt_vec_vis*pt_vec_invis)))
    return mT 

# Sample Inputs 
bottom_sideA_array = np.array([1, 2, 3, 4]) # (px, py, pz, E)
lepton_sideA_array = np.array([1, 1, 2, 5]) 
vis_sideA_array = np.array([bottom_sideA_array, lepton_sideA_array]) 

bottom_sideB_array = np.array([1, 1, 1, 4]) # (px, py, pz, E)
lepton_sideB_array = np.array([1, 3, 4, 5]) 
vis_sideB_array = np.array([bottom_sideB_array, lepton_sideB_array]) 

met = np.array([1, 3, 0, 5]) # (px, py, 0, E)
alphaList = [0.9] 

# Function to Minimise 
def objective(invis_sideA_array): 
    """The function to minimise is the linear combination of alpha*m(vis_parent_layer1)' + beta*m(vis_parent_layer2)'."""
    alpha_term_1 = mT_4vecCalc(vis_sideA_array[-1], invis_sideA_array) # mT(lA, pT_A)
    alpha_term_2 = mT_4vecCalc(vis_sideB_array[-1], met-invis_sideA_array) # mT(TB, pT_B) 
    alpha_term = max(alpha_term_1, alpha_term_2) 
    
    beta_term_1 = mT_4vecCalc(vis_sideA_array[0] + vis_sideA_array[-1], invis_sideA_array) # mT(lATA, pT_A)
    beta_term_2 = mT_4vecCalc(vis_sideB_array[0] + vis_sideB_array[-1], met-invis_sideA_array) # mT(TBbB, pt_B)
    beta_term = max(beta_term_1, beta_term_2) 
    
    return alphaList[0]*alpha_term + (1-alphaList[0])*beta_term 

# Example guess & minimisation execution 
invis_side_A_array_guess = np.array([1, 2, 0, 1]) # [px, py, pz, E] guess of the invisible neutrino of side A 

#bnds = [(None, None), (None, None), (0, 0), (None, None)] # bounds 
sol = so.minimize(objective, x0 = invis_side_A_array_guess, method='SLSQP') 

# TESTS OF FUNCTIONALITY 
# Et_scalarCalc(m, pt) 
Et_scalarCalc(3, 4) == 5 # True 

# mT_4vecCalc(p4_vis_array, p4_invis_array)
p4_vis_array = np.array([1, 2, 3, 4])
p4_vis = ROOT.TLorentzVector() 
p4_vis.SetPxPyPzE(p4_vis_array[0], p4_vis_array[1], p4_vis_array[2], p4_vis_array[3])

p4_invis_array = np.array([1, 2, 4, 6])
p4_invis = ROOT.TLorentzVector() 
p4_invis.SetPxPyPzE(p4_invis_array[0], p4_invis_array[1], p4_invis_array[2], p4_invis_array[3])

pt_vec_vis = p4_vis.Vect() 
pt_vec_vis.SetZ(0)
pt_vec_invis = p4_invis.Vect()
Et_vis = Et_scalarCalc(p4_vis.M(), p4_vis.Pt()) 
Et_invis = Et_scalarCalc(p4_invis.M(), p4_invis.Pt()) 

p4_vis.Px() == 1 # True 
p4_vis.Py() == 2 # True 
p4_vis.Pz() == 3 # True 
p4_vis.E() == 4 # True
p4_vis.M() == math.sqrt(4**2 - (1+2**2 + 3**2)) # True 
p4_vis.Pt() == math.sqrt(1**2 + 2**2) # True 
p4_vis.Phi() == 1.1071487177940904 # ~np.arctan(2), True 
Et_vis == 2.6457513110645907 # True 
pt_vec_vis.Px() == 1 # True 
pt_vec_vis.Py() == 2 # True
pt_vec_vis.Pz() == 0 # True 
pt_vec_vis.Pt() == math.sqrt(1**2 + 2**2) # True 
mT_4vecCalc(p4_vis_array, p4_invis_array) == 5.537537280452248 # True 

# optimise
# will return non-zero pz value for solution vector, negative px, py value and v. small E value (not probable) ?? 
