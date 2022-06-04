###############################################
# Mt2dc (mt2 decay chain)
# Purpose: Calculate the kinematic variable mt2dc.
#    This requires doing a numerical minimization.
#    The minimisation determines the pt of two
#    dummy vectors that represent the invisible 
#    particles in the event.

# Event considered (slide 9) 
############################################### 

import ROOT
import numpy as np
import math 
from scipy.optimize import minimize

def Et_scalarCalc(m, pt):
    """Calculate the transverse energy of a particle, using scalar input variables."""
    Et = math.sqrt(max(0, m**2+pt**2)) 
    return Et

def mT_scalarcalc(m_vis, m_invis, pt_vis, pt_invis, cosAngle):
    """Calculate the transverse mass of a parent particle, using scalar input variables.
    In this case, a parent particle decays to produce visible and invisible products."""
    Et_vis = Et_scalarCalc(m_vis, pt_vis)
    Et_invis = Et_scalarCalc(m_invis, pt_invis) 
    
    mT_squared = max(0, m_vis*m_vis + m_invis*m_invis + 2*(Et_vis*Et_invis - pt_vis*pt_invis*cosAngle) )
    return math.sqrt(mT_squared)
           
    
def mT_4vecCalc(p4_vis, p4_invis):
    """Calculate the tranverse mass of a parent particle, using 4-vector input variables. 
    In this case, a parent particle decays to produce visible and invisible products. 
    """
    m_vis = max(0, p4_vis.M()) # get M
    m_invis = max(0, p4_invis.M()) 
    Et_vis = Et_scalarCalc(m_vis, p4_vis.Pt()) 
    Et_invis = Et_scalarCalc(m_invis, p4_invis.Pt()) 
    
    # Make spatial 3-vectors for momentum calculation 
    pt_vec_vis = p4_vis.Vect() # (E, px, py, pz) -> (px, py, pz)
    pt_vec_vis.SetZ(0) # then, set pz = 0, to obtain a transverse momemtnum vector 
    pt_vec_invis = p4_invis.Vect() 
    pt_vec_invis.SetZ(0) 
    
    mT = math.sqrt(max(0, m_vis**2 + m_invis**2 + 2*(Et_vis*Et_invis - pt_vec_vis*pt_vec_invis)))
    return mT 
       
class mt2dc: 
   # Defining the types of relevant class attributes 
   TLV_temp = ROOT.TLorentzVector() # (E, px, py, pz) 
   _TLV_vis_sideA = np.array([TLV_temp]) 
   _TLV_vis_sideB = np.array([TLV_temp]) # array to store all TLV of vis. particles on side A 
    # order: left layer (e.g., bottom A), right layer (e.g., lepton A
    # order: left layer (e.g., bottom B), right layer (e.g., lepton B)
   _met = ROOT.TLorentzVector() # missing transverse E (in which is comprehended missing transverse pt?) 
   _mass_inv_sideA = 0
   _mass_inv_sideB = 0
   _alphaList = np.array([0]) # array to contain the list of scale factors,  e.g. [0.7, 0.3]

   # Required in the minimisation 
   _gamma_sideA_2vec = ROOT.TVector2() # (px, py) of the substitute particle (for neutrino) of side A 
   _gamma_sideB_2vec = ROOT.TVector2() # (px, py) of the substitute particle (for neutrino) of side B 
    
   # Resulting quantities 
   _mt2prime_list = np.array([1]) # fill with individual mt2' (e.g., 2 elements if 2-layer study)
   _mT2dc = -1 # final value acquired 
   
   _calcPerformed = False

   def __init__(self, alphaList, met, TLV_List_vis_sideA, TLV_List_vis_sideB, gamma_sideA_2vec, 
                gamma_sideB_2vec, mass_inv_sideA, mass_inv_sideB, mt2prime_list):
      self._alphaList = np.array(alphaList)
      self._met = met 
      self._TLV_List_vis_sideA = TLV_List_vis_sideA
      self._TLV_List_vis_sideB = TLV_List_vis_sideB
      self._gamma_sideA_2vec = gamma_sideA_2vec 
      self._gamma_sideB_2vec = gamma_sideB_2vec 
      self._mass_inv_sideA = mass_inv_sideA
      self._mass_inv_sideB = mass_inv_sideB 
      self._mt2prime_list = mt2prime_list 
      
      if (  ( self._alphaList.size != self._TLV_vis_sideA.size ) or   \
            ( self._alphaList.size != self._TLV_vis_sideB.size )          ):
         print("List lengths are not all the same")
         print("self._alphaList.size = ", self._alphaList.size)
         print("self._V_layersList_sideA.size = ", self._TLV_vis_sideA.size)
         print("self._V_layersList_sideB.size = ", self._TLV_vis_sideB.size)
      pass
      
   def get_alphaList(self):
      return self._alphaList.tolist()
   def get_met(self):
      return self._met
   def get_TLV_List_vis_sideA(self):
      return self._TLV_List_vis_sideA 
   def get_TLV_List_vis_sideB(self):
      return self._TLV__Listvis_sideB
   def get_DummyInvis_A_2vec(self):
      return self._gamma_sideA_2vec
   def get_DummyInvis_B_2vec(self):
      return self._gamma_sideB_2vec
   def get_alpha_mT2_prime_layersList(self):
      return (self._alphaList * self._mt2prime_list).tolist()
    # multiply each alpha by the mt2_prime value to which it corresponds; used to convert to a list 
   def get_mass_inv_sideA(self):
      return self._mass_inv_sideA 
   def get_mass_inv_sideB(self):
      return self._mass_inv_sideB 
   def get_mT2_prime_layersList(self):
      return self._mt2prime_list.tolist()
   def get_mT2dc(self):
      if not self._calcPerformed:
         print("WARNING: Minimization not yet performed!")
         return -1
      return self._mT2dc
   def get_muT2dc(self, mass_parents_list):
      """Get the normalised muT2DC quantity."""
      mMax = self.get_mT2dc_maxMass(mass_parents_list)
      return self.get_mT2dc() / mMax

   def get_mT2dc_maxMass(self, mass_parents_list):
       """Get the upper bound expected of the mT2dc mass calculation."""
       mParentsList = np.array(mass_parents_list)
       if (  self._alphaList.size != mParentsList.size ):
          print("List lengths are not all the same")
          print("self._alphaList.size = ", self._alphaList.size)
          print("length of mass_parents_list = ", mParentsList.size)
          return -1
       return sum(self._alphaList * mParentsList)

   def get_allSettings(self):
      """Function that returns all settings inputted by user, required for mT2DC
      calculation. Quantities are all in a single list so that user can make sure
      they all look correct at the same time with less risk of calling the
      information being incorrectly recorded from calling each of the
      individual get functions."""
    
      allSettingsList = [  self._TLV_List_vis_sideA.tolist(),    \
                           self._TLV_List_vis_sideB.tolist(),    \
                           self._met,                            \
                           self._mass_inv_sideA, self._mass_inv_sideB,   \
                           self._alphaList.tolist()                  ]
      return allSettingsList

   def get_allSettingsAndResults(self):
        """Function that returns all the settings inputted by the user for
        doing the mt2dc calculation as well as the results."""
        allSettingsResultsList =  [self.get_allSettings(),
                                  [ self._mT2dc,
                                     self._mt2prime_list.tolist(),
                                     self._gamma_sideA_2vec,
                                     self._gamma_sideB_2vec,
                                     self._calcPerformed ]]
        return allSettingsResultsList
    
   def calculate_mt2dc(self):
        """Function to execute the calculation of mt2dc.
        """
        # for a two-layer study, observe that 
        # _TLV_vis_sideA = [vis_leftmost_layer – b, vis_rightmost_layer – l]
        # alpha_term = mT_4vecCalc(p4_vis - l, p4_invis - neutrino)
        
        def objective(TLV_invis_sideA, TLV_invis_sideB): 
            """Defined exterior to this function:  TLV_List_vis_sideA, TLV_invis_sideB, alpha.
            TLV_invis_sideA = [E_T_sideA, px_A, py_A, 0] 
            TLV_invis_sideB = [E_T_sideB, px_B, py_B, 0] 
            """
            
            alpha_term_1 = mT_4vecCalc(_TLV_vis_sideA[-1], TLV_invis_sideA) # mT(lA, pT_A)
            alpha_term_2 = mT_4vecCalc(_TLV_vis_sideB[-1], TLV_invis_sideB) # mT(TB, pT_B) 
            alpha_term = max(alpha_term_1, alpha_term_2) 
            
            beta_term_1 = mT_4vecCalc(_TLV_vis_sideA[0] + _TLV_vis_sideA[-1], TLV_invis_sideA) # mT(lATA, pT_A)
            beta_term_2 = mT_4vecCalc(_TLV_vis_sideB[0] + _TLV_vis_sideB[-1], TLV_invis_sideB) # mT(TBbB, pt_B)
            beta_term = max(beta_term_1, beta_term_2) 
            
            return alphaList[0]*alpha_term + (1-alphaList[0])*beta_term 
         
       def constraint(TLV_invis_sideA, TLV_invis_sideB): 
            """pT_A + pT_B - pT_miss = 0.
            Supposing met = [ET_miss, pT_miss, 0] Lorentz 4-vector
            """
            
            pT_A = TLV_invis_sideA.Pt()
            pT_B = TLV_invis_sideB.Pt() 
            
            return pT_A - pT_B - met.Pt() 
        
       con1 = {'type': 'eq', 'fun': constraint}
    
       sol = minimize(objective, constraints=cons) 
        
        _calcPerformed = True
        if True:
            return True
        else:
            return False
    
    
