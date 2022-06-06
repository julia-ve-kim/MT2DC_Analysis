###############################################
# Mt2dc: Calculate the kinematic variable mt2dc.
# This requires doing a numerical minimization.
# The minimisation determines the pt of two
# dummy vectors that represent the invisible 
# particles in the event.
###############################################

import ROOT
import numpy as np 
import math 
import scipy.optimize as so 

# FUNCTIONS 
# Basic Variable Functions 
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

# TLorentz Module Functions 
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

# Optimisation Functions 
def get_alpha_term(vis_sideA_array, vis_sideB_array, met, invis_sideA_array_variable):
    alpha_term_1 = mT_arrayCalc(vis_sideA_array[-1], invis_sideA_array_variable) # mT(lA, pT_A)
    alpha_term_2 = mT_arrayCalc(vis_sideB_array[-1], met-invis_sideA_array_variable) # mT(TB, pT_B) 
    alpha_term = max(alpha_term_1, alpha_term_2)
    return alpha_term 