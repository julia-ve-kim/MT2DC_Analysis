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
       
class mt2dc: 
    # initialise types 
    _TLV_array_temp = np.array([])
    _vis_sideA_array = np.array([_TLV_array_temp]) 
    _vis_sideB_array = np.array([_TLV_array_temp]) 
    _invis_sideA_array = np.array([]) 
    _invis_sideB_array = np.array([]) 
    _met = np.array([_TLV_array_temp]) 
    _alphaList = [] 
    _gamma_sideA_2vec = ROOT.TVector2() # (px, py) 
    _gamma_sideB_2vec = ROOT.TVector2() 
    _mass_inv_sideA = 0
    _mass_inv_sideB = 0 
    _invis_side_A_array_guess = np.array([]) 
    _calcPerformed = False 
    
    # resulting quantities 
    _mt2prime_list = [] # filled with mt2' 
    _mt2dc = 0 # final value acquired 
    
    def __init__(self, vis_sideA_array, vis_sideB_array, invis_sideA_array, invis_sideB_array, met, alphaList, gamma_sideA_2vec, 
                 gamma_sideB_2vec, mass_inv_sideA, mass_inv_sideB, invis_side_A_array_guess, mt2prime_list, mt2dc): 
        self._vis_sideA_array = vis_sideA_array
        self_vis_sideB_array = vis_sideB_array 
        self._invis_sideA_array = invis_sideA_array 
        self._invis_sideB_array = invis_sideB_array 
        self._met = met 
        self._alphaList = alphaList 
        self._gamma_sideA_2vec = gamma_sideA_2vec 
        self._gamma_sideB_2vec = gamma_sideB_2vec 
        self._mass_inv_sideA = mass_inv_sideA
        self._mass_inv_sideB = mass_inv_sideB 
        self._invis_side_A_array_guess = invis_side_A_array_guess 
        self._mt2prime_list = mt2prime_list
        self._mt2dc = mt2dc 
      
        if ((self._alphaList.size != self._vis_sideA_array.size) or  \
            (self._alphaList.size != self._vis_sideB_array.size)):
            print("List lengths are not all the same")
            print("self._alphaList.size = ", self._alphaList.size)
            print("self._vis_sideA_array.size = ", self._vis_sideA_array.size)
            print("self._vis_sideB_array.size = ", self._vis_sideB_array.size)
            pass
        
        # function attributes 
        def get_vis_sideA_array(self):
            return self._vis_sideA_array 
        
        def get_vis_sideB_array(self):
            return self._vis_sideB_array
        
        def get_alphaList(self):
            return self._alphaList.tolist() 
        
        def get_met(self):
            return self._met 

        def get_mass_inv_sideA(self):
            return self._mass_inv_sideA 
        
        def get_mass_inv_sideB(self):
            return self._mass_inv_sideB 
        
        def get_invis_side_A_array_guess(self):
            return self._invis_side_A_array_guess  
        
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
            alpha_term_2 = mT_4vecCalc(vis_sideB_array[-1], met-invis_sideA_array_variable) # mT(TB, pT_B) 
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
