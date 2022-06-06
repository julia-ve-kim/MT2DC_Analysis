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
from scipy.optimize import minimize 

# test values
# vis_sideA_array = np.array([[1, 2, 3, 4], [1, 2, 3, 4]]) 
# vis_sideB_array = np.array([[1, 2, 3, 4], [1, 2, 3, 4]]) 
# invis_sideA_array = np.array([1, 2, 3, 4])
# invis_sideB_array = np.array([1, 2, 3, 4]) 
# met = np.array([1, 2, 3, 4]) 
# alphaList = np.array([1]) 
# gamma_sideA_2vec = np.array([1, 2]) 
# gamma_sideB_2vec = np.array([3, 4]) 
# mass_inv_sideA = 1
# mass_inv_sideB = 2 
# invis_side_A_array_guess = np.array([1, 2, 3, 4]) 
# mt2prime_list = np.array([1, 2])
# mt2dc = 1 

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
 
class mt2dc:
    """This class presumes the minimisation function has already been run...? 
    """
    _calcPerformed = False 
    
    def __init__(self, _vis_sideA_array, _vis_sideB_array, _invis_sideA_array, _invis_sideB_array, _met, _alphaList,
                 _gamma_sideA_2vec, _gamma_sideB_2vec, _mass_inv_sideA, _mass_inv_sideB, _invis_side_A_array_guess, 
                 _mt2prime_list, _mt2dc_val): 
        # the following enable you to "get" the quantities, by using the syntax a = mt2dc(...), a._quantity 
        self._vis_sideA_array = _vis_sideA_array
        self._vis_sideB_array = _vis_sideB_array 
        self._invis_sideA_array = _invis_sideA_array 
        self._invis_sideB_array = _invis_sideB_array 
        self._met = _met 
        self._alphaList = _alphaList 
        self._gamma_sideA_2vec = _gamma_sideA_2vec 
        self._gamma_sideB_2vec = _gamma_sideB_2vec 
        self._mass_inv_sideA = _mass_inv_sideA
        self._mass_inv_sideB = _mass_inv_sideB 
        self._invis_side_A_array_guess = _invis_side_A_array_guess 
        self._mt2prime_list = _mt2prime_list
        self._mt2dc = _mt2dc_val
      
        if ((self._alphaList.size != self._vis_sideA_array.size) or  \
            (self._alphaList.size != self._vis_sideB_array.size)):
            print("List lengths are not all the same")
            print("self._alphaList.size = ", self._alphaList.size)
            print("self._vis_sideA_array.size = ", self._vis_sideA_array.size)
            print("self._vis_sideB_array.size = ", self._vis_sideB_array.size)
            pass
        
        # function attributes  
        def run_objective(self, invis_sideA_array_variable): 
            """This function will serve as the objective function, being equal to the inner expression of mt2dc."""
            alpha_term_1 = mT_4vecCalc(vis_sideA_array[-1], invis_sideA_array_variable) # mT(lA, pT_A)
            alpha_term_2 = mT_4vecCalc(vis_sideB_array[-1], met-invis_sideA_array_variable) # mT(TB, pT_B) 
            alpha_term = max(alpha_term_1, alpha_term_2) 
            
            beta_term_1 = mT_4vecCalc(vis_sideA_array[0] + vis_sideA_array[-1], invis_sideA_array_variable) # mT(lATA, pT_A)
            beta_term_2 = mT_4vecCalc(vis_sideB_array[0] + vis_sideB_array[-1], met-invis_sideA_array_variable) # mT(TBbB, pt_B)
            beta_term = max(beta_term_1, beta_term_2) 
            
            return alphaList[0]*alpha_term + (1-alphaList[0])*beta_term 
        
        def get_alpha_term(self, invis_sideA_array_variable):
            alpha_term_1 = mT_4vecCalc(vis_sideA_array[-1], invis_sideA_array_variable) # mT(lA, pT_A)
            alpha_term_2 = mT_4vecCalc(vis_yousideB_array[-1], met-invis_sideA_array_variable) # mT(TB, pT_B) 
            alpha_term = max(alpha_term_1, alpha_term_2)
            return alpha_term 
        
        def get_beta_term(self, invis_sideA_array_variable):
            beta_term_1 = mT_4vecCalc(vis_sideA_array[0] + vis_sideA_array[-1], invis_sideA_array_variable) # mT(lATA, pT_A)
            beta_term_2 = mT_4vecCalc(vis_sideB_array[0] + vis_sideB_array[-1], met-invis_sideA_array_variable) # mT(TBbB, pt_B)
            beta_term = max(beta_term_1, beta_term_2) 
            return beta_term 
            
        def calculate_mt2dc(self, invis_side_A_array_guess):
            """Function to execute the calculation of mt2dc."""
            sol = so.minimize(run_objective, x0 = invis_side_A_array_guess, method='SLSQP') # [px, py, pz, E]
            
            if sol.success: 
                self._calcPerformed = True 
                return sol 
            else:
                return -1 
         
        def get_invis_sideA_array(self):
            sol = calculate_mt2dc(self, invis_side_A_array_guess)
            self._invis_sideA_array = sol.x 
            return self._invis_sideA_array 
        
        def get_invis_sideB_array(self):
            sol = calculate_mt2dc(self, invis_side_A_array_guess)
            if self._calcPerformed:
                self._invis_sideB_array = met - sol.x 
                return self._invis_sideB_array
            else:
                print("Minimisation has not been successfully performed!") 
                return -1  
       
        def get_DummyInvis_A_2vec(self):
            sol = calculate_mt2dc(self, invis_side_A_array_guess)
            if self._calcPerformed:
                self._gamma_sideA_2vec = ROOT.TVector2(sol.x[0], sol.x[1]) 
                return self._gamma_sideA_2vec 
            else:
                print("Minimisation has not been successfully performed!") 
                return -1 
        
        def get_DummyInvis_B_2vec(self):
            sol = calculate_mt2dc(self, invis_side_A_array_guess)
            if self._calcPerformed: 
                self._gamma_sideB_2vec = ROOT.TVector2(met[0]-sol.x[0], met[1]-sol.x[1]) 
                return self._gamma_sideB_2vec
            else:
                print("Minimisation has not been successfully performed!")
                return -1 
        
        def get_mt2_prime_layersList(self): 
            mt2prime_list.append(get_alpha_term(self, invis_side_A_array))
            mt2prime_list.append(get_beta_term(self, invis_side_A_array))
            
            self._mt2prime_list = mt2prime_list 
            return self._mt2prime_list.tolist()
        
        def get_alpha_mT2_prime_layersList(self):
            return (self._alphaList*self._mt2prime_list).tolist()
        
        def get_mt2dc(self): 
            sol = calculate_mt2dc(self, invis_side_A_array_guess)
            if self._calcPerformed: 
                self._mt2dc = sol.fun 
                return self._mt2dc
            else:
                print("Minimisation has not been successfully performed!")
                return -1 
        
        def get_mt2dc_maxMass(self, mass_parents_list):
            """Get the upper bound expected of the mT2dc mass calculation."""
            mParentsList = np.array(mass_parents_list)
            if (self._alphaList.size != mParentsList.size):
                print("List lengths are not all the same")
                print("self._alphaList.size = ", self._alphaList.size)
                print("length of mass_parents_list = ", mParentsList.size)
                return -1
            return sum(self._alphaList * mParentsList)
        
        def get_mut2dc(self, mass_parents_list):
            """Get the normalised mut2DC quantity."""
            mMax = self.get_mt2dc_maxMass(mass_parents_list)
            return self.get_mt2dc() / mMax
        
        def get_allSettings(self):
            """Function that returns all settings inputted by user, required for mT2DC
            calculation. Quantities are all in a single list so that user can make sure
            they all look correct at the same time with less risk of calling the
            information being incorrectly recorded from calling each of the
            individual get functions."""
            allSettingsList = [self._vis_sideA_array.tolist(),    \
                           self._vis_sideB_array.tolist(),    \
                           self._met,                            \
                           self._mass_inv_sideA, self._mass_inv_sideB,   \
                           self._alphaList.tolist()]
            return allSettingsList
        
        def get_allSettingsAndResults(self):
            """Function that returns all the settings inputted by the user for
            doing the mt2dc calculation as well as the results."""
            allSettingsResultsList =  [self.get_allSettings(),
                                  [self._mt2dc,
                                     self._mt2prime_list.tolist(),
                                     self._gamma_sideA_2vec,
                                     self._gamma_sideB_2vec,
                                     self._calcPerformed]]
            return allSettingsResultsList
