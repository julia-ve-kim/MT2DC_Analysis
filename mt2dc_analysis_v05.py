###############################################
# Mt2dc Analysis
##############################################

import ROOT
import mt2dc_v02_prime as DC 
import math
import numpy as np 
import scipy.optimize as so 
from array import array 

##############################################
# Define input & output files 
##############################################
f_inputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_03_March_07_skim_mg5_ttbar_jet_merged_001-716_ntuple_2l2b_v01.root", "read")
t = f_inputRoot.Get("variables")
type(t)
outDir = "/Users/juliakim/Documents/mT2DCAnalysisPlots/"  

##############################################
# Define tree 
##############################################
f_outputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_05_May_10_mt2dc_analysis_v10.root", "recreate")
tree = ROOT.TTree("results", "tree storing mt2dc calculation results") 

# create float array pointers 
mT2dc_diff = array('d', [0]) # sign(+-) => may be used in weighted averages calculation 
mT2dc = array('d', [0])
mT2_W = array('d', [0]) 
mT2prime_W = array('d', [0]) 
mT2prime_W_subC = array('d', [0]) # consider only mT2_W < 1 events 
mT2_t = array('d', [0]) 
mT2prime_t = array('d', [0]) 
mT2prime_t_subC = array('d', [0]) # consider only mT2_W < 1 events 
sub_pT_sideA = array('d', [0]) 
sub_pT_sideB = array('d', [0]) 
sub_pT_min_over_met = array('d', [0]) # min(sub_pT_sideA, sub_pT_sideB)/met
alpha = array('d', [0])
constraint_pT_cut = array('d', [0]) 
constraint_pT_subcut[0] = array('d', [0]) 
calc_speed = array('d', [0]) # 1 = fast, 0 = slow 
success = array('d', [0]) # 1 = success, 0 = failure 

# create branches
tree.Branch("mT2dc_diff", mT2dc_diff, 'mT2dc_diff/D') 
tree.Branch("mT2dc", mT2dc, 'mT2dc/D') 
tree.Branch("mT2_W", mT2_W, 'mT2_W/D') 
tree.Branch("mT2prime_W", mT2prime_W, 'mT2prime_W/D') 
tree.Branch("mT2prime_W_subC", mT2prime_W_subC, 'mT2prime_W_subC/D') 
tree.Branch("mT2_t", mT2_t, 'mT2_t/D') 
tree.Branch("mT2prime_t", mT2prime_t, 'mT2prime_t/D') 
tree.Branch("mT2prime_t_subC", mT2prime_t_subC, 'mT2prime_t_subC/D') 
tree.Branch("sub_pT_sideA", sub_pT_sideA, 'sub_pT_sideA/D') 
tree.Branch("sub_pT_sideB", sub_pT_sideB, 'sub_pT_sideB/D') 
tree.Branch("sub_pT_min_over_met", sub_pT_min_over_met, 'sub_pT_min_over_met/D') 
tree.Branch("alpha", alpha, 'alpha/D') 
tree.Branch("constraint_pT_cut", constraint_pT_cut, 'constraint_pT_cut/D') 
tree.Branch("calc_speed", calc_speed, 'calc_speed/D') 
tree.Branch("success", success, 'success/D') 
                
##############################################
# Create plots of variables in input ROOT file 
##############################################         
h_ell1_pt = ROOT.TH1F("h_ell1_pt", "Pt of highest pt lepton; Leading light lepton pt [GeV]; #entries/2 GeV", 100, 0, 200)
h_ell1_E = ROOT.TH1F("h_ell1_E", "E of highest pt lepton; Leading light lepton Energy [GeV]; #entries/2.5 GeV", 100, 0, 250)

h_ell2_pt = ROOT.TH1F("h_ell2_pt", "Pt of lowest pt lepton; Second light lepton pt [GeV]; #entries/2GeV", 100, 0, 200)
h_ell2_E = ROOT.TH1F("h_ell2_pt", "E of lowest pt lepton; Second light lepton pt [GeV]; #entries/2.5GeV", 100, 0, 250)

h_bjet1_E = ROOT.TH1F("h_bjet1_E", "E of highest pt b-tagged jet; Leading b-tagged jet Energy [GeV]; #entries/5GeV", 100, 0, 500)
h_bjet2_E = ROOT.TH1F("h_bjet_2_E", "E of lowest pt b-tagged jet; Leading b-tagged jet Energy [GeV]; #entries/5GeV", 100, 0, 500)

h_mT2_W = ROOT.TH1F("h_mT2_W", "mt2(ell1,ell2) = mt2(W); mt2(W) [GeV]; #entries/1GeV", 200, 0, 200)
h_mT2_t_11_22 = ROOT.TH1F("h_mT2_t_11_22", "mt2|11,22(t) = mt2(b1 ell1,b2 ell2); mt2(t|11,22) [GeV]; #entries/1GeV", 300, 0, 300)
h_mT2_t_12_21 = ROOT.TH1F("h_mT2_t_12_21", "mt2|12,21(t) = mt2(b1 ell2,b2 ell1); mt2(t|12,21) [GeV]; #entries/1GeV", 300, 0, 300)
h_mT2_t_min = ROOT.TH1F("h_mT2_t_min", "min( mt2(t|11,22), mt2(t|12,21)); mt2(t)_min [GeV]; #entries/1GeV", 300, 0, 300)

h_EtMiss = ROOT.TH1F("h_EtMiss", "Etmiss; Etmiss [GeV]; #entries/1GeV", 300, 0, 300)
h_EtMiss_phi = ROOT.TH1F("h_EtMiss_phi", "Etmiss azimuthal angle; Azimuthal angle [rad]; #entries/0.02rad", 400, -4, 4)  

##############################################
# Define constants and optimisation mode 
##############################################
m_W = 80 # GeV 
m_t = 173 # GeV
nentries = t.GetEntries() 
calcStyle = input("Enter 'fast' or 'slow':") 
alphaList = np.linspace(0, 1, 20)  
pT_cut = 20 # GeV 

##############################################
# Main analysis - loop over all events
##############################################
for i in range(20000):
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
    
    # RETRIEVE INFORMATION FROM TREE TO FILL HISTOGRAMS. 
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
    
    # CALCULATION OF MT2DC:
    # side A neutrino pT vector guesses 
    invis_sideA_array_guess_1 = met[:2]/2 
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
    
    invis_sideA_array_guesses = [invis_sideA_array_guess_1, invis_sideA_array_guess_2, invis_sideA_array_guess_3, 
                                 invis_sideA_array_guess_4, invis_sideA_array_guess_5, invis_sideA_array_guess_6, 
                                 invis_sideA_array_guess_7, invis_sideA_array_guess_8, invis_sideA_array_guess_9, 
                                 invis_sideA_array_guess_10, invis_sideA_array_guess_11, invis_sideA_array_guess_12,
                                 invis_sideA_array_guess_13, invis_sideA_array_guess_14, invis_sideA_array_guess_15,
                                 invis_sideA_array_guess_16, invis_sideA_array_guess_17, invis_sideA_array_guess_18, 
                                 invis_sideA_array_guess_19, invis_sideA_array_guess_20, invis_sideA_array_guess_21]
    
    def constraint_1(invis_sideA_2vec):
        return invis_sideA_2vec[0]**2 + invis_sideA_2vec[1]**2 - pT_cut**2 
    
    def constraint_2(invis_sideA_2vec):
        invis_sideB_2vec = met[:2] - invis_sideA_2vec 
        return invis_sideB_2vec[0]**2 + invis_sideB_2vec[1]**2 - pT_cut**2 
    
    cons = [{'type': 'ineq', 'fun': constraint_1},  {'type': 'ineq', 'fun': constraint_2}] 
        
    if calcStyle == 'fast':
        for alphaVal in alphaList: 
            def objective(invis_sideA_2vec): 
                invis_sideA_array = np.array([invis_sideA_2vec[0], invis_sideA_2vec[1], 0, 
                                              np.sqrt(invis_sideA_2vec[0]**2 + invis_sideA_2vec[0]**2)]) 
                
                alpha_term_1 = DC.mT_arrayCalc(vis_sideA_array[-1], invis_sideA_array) 
                alpha_term_2 = DC.mT_arrayCalc(vis_sideB_array[-1], met-invis_sideA_array) 
                alpha_term = max(alpha_term_1, alpha_term_2) 
                
                beta_term_1 = DC.mT_arrayCalc(vis_sideA_array[0] + vis_sideA_array[-1], invis_sideA_array) 
                beta_term_2 = DC.mT_arrayCalc(vis_sideB_array[0] + vis_sideB_array[-1], met-invis_sideA_array) 
                beta_term = max(beta_term_1, beta_term_2) 
                
                return alphaVal*alpha_term + (1-alphaVal)*beta_term 
            
            # evaluate objective function at each guess 
            objective_at_guess = [] 
            for guess in invis_sideA_array_guesses: 
                objective_at_guess.append(objective(guess)) 
            
            # optimise objective function at guess for which objective_at_guess is smallest 
            sol_UC = so.minimize(objective, x0 = invis_sideA_array_guesses[np.argmin(objective_at_guess)], method='SLSQP', 
                                 options={'maxiter': 2000, 'ftol': 1e-7,'disp': True}) # unconstrained opt. 
            sol = so.minimize(objective, x0 = invis_sideA_array_guesses[np.argmin(objective_at_guess)], method='SLSQP', 
                              options={'maxiter': 2000, 'ftol': 1e-7,'disp': True}, constraints=cons) # constraint opt. 
            
            # > CONSTRAINED 
            if sol.success == True: 
                # filling tree  
                mT2dc_diff[0] = sol.fun - (alphaVal*mt2_W + (1-alphaVal)*mt2_t_11_22)
                mT2dc[0] = sol.fun
                mT2prime_W[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol.x)
                mT2prime_t[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol.x)
                sub_pT_sideA[0] = np.linalg.norm(sol.x)
                sub_pT_sideB[0] = np.linalg.norm(met[:2] - sol.x)
                sub_pT_min_over_met[0] = (min(np.linalg.norm(sol.x), np.linalg.norm(met[:2] - sol.x))/np.linalg.norm(met[:2]))
                alpha[0] = alphaVal 
                constraint_pT_cut[0] = pT_cut  
                calc_speed[0] = 1 
                success[0] = 1
                
                tree.Fill() 
                
                if mT2_W < 1:
                    mT2prime_W_subC[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol.x) 
                    mT2prime_t_subC[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol.x) 
                    constraint_pT_subcut[0] = pT_cut 
                    tree.Fill() 
                
            else:
                # optimise objective function using COBYLA in lieu of SLSQP (which encounters no exit mode 8 error) 
                sol_alt = so.minimize(objective, x0 = invis_sideA_array_guesses[np.argmin(objective_at_guess)], method='COBYLA', 
                          options={'maxiter': 2000, 'ftol': 1e-7,'disp': True}, constraints=cons) # alternative solution  
                # filling tree 
                mT2dc_diff[0] = sol_alt.fun - (alphaVal*mt2_W + (1-alphaVal)*mt2_t_11_22)
                mT2dc[0] = sol_alt.fun
                mT2prime_W[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_alt.x)
                mT2prime_t[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_alt.x)
                sub_pT_sideA[0] = np.linalg.norm(sol_alt.x)
                sub_pT_sideB[0] = np.linalg.norm(met[:2] - sol_alt.x)
                sub_pT_min_over_met[0] = (min(np.linalg.norm(sol_alt.x), np.linalg.norm(met[:2] - 
                                                                                        sol_alt.x))/np.linalg.norm(met[:2]))
                alpha[0] = alphaVal 
                constraint_pT_cut[0] = pT_cut  
                calc_speed[0] = 1 
                
                if mT2_W < 1:
                    mT2prime_W_subC[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_alt.x) 
                    mT2prime_t_subC[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_alt.x) 
                    constraint_pT_subcut[0] = pT_cut 
                    tree.Fill() 
                
                if sol_alt.success == True:
                    success[0] = 1
                else:
                    success[0] = 0 
                
                tree.Fill() 
            
            # > UNCONSTRAINED 
            mT2dc_diff[0] = sol_UC.fun - (alphaVal*mt2_W + (1-alphaVal)*mt2_t_11_22)
            mT2dc[0] = sol_UC.fun
            mT2prime_W[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x)
            mT2prime_t[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x)
            sub_pT_sideA[0] = np.linalg.norm(sol_UC.x)
            sub_pT_sideB[0] = np.linalg.norm(met[:2] - sol_UC.x)
            sub_pT_min_over_met[0] = (min(np.linalg.norm(sol_UC.x), np.linalg.norm(met[:2] - 
                                                                                        sol_UC.x))/np.linalg.norm(met[:2]))
            mT2_W[0] = mt2_W 
            mT2_t[0] = mt2_t 
            alpha[0] = alphaVal 
            constraint_pT_cut[0] = 0 
            calc_speed[0] = 1 
            success[0] = 1
            
            if mT2_W < 1:
                mT2prime_W_subC[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x) 
                mT2prime_t_subC[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x) 
                constraint_pT_subcut[0] = 0 
                tree.Fill() 
                         
    elif calcStyle == 'slow':
        for alphaVal in alphaList: 
            def objective(invis_sideA_2vec): 
                invis_sideA_array = np.array([invis_sideA_2vec[0], invis_sideA_2vec[1], 0, 
                                              np.sqrt(invis_sideA_2vec[0]**2 + invis_sideA_2vec[0]**2)]) 
                
                alpha_term_1 = DC.mT_arrayCalc(vis_sideA_array[-1], invis_sideA_array) 
                alpha_term_2 = DC.mT_arrayCalc(vis_sideB_array[-1], met-invis_sideA_array) 
                alpha_term = max(alpha_term_1, alpha_term_2) 
                
                beta_term_1 = DC.mT_arrayCalc(vis_sideA_array[0] + vis_sideA_array[-1], invis_sideA_array) 
                beta_term_2 = DC.mT_arrayCalc(vis_sideB_array[0] + vis_sideB_array[-1], met-invis_sideA_array) 
                beta_term = max(beta_term_1, beta_term_2) 
                
                return alphaVal*alpha_term + (1-alphaVal)*beta_term 
            
            # run minimiser using every initial guess in invis_side_array_guess
            minimised_objective_at_guess_UC = [] 
            minimised_objective_at_guess_UC_fun = [] 
            minimised_objective_at_guess = [] 
            minimised_objective_at_guess_fun = [] 
            
            for guess in invis_sideA_array_guesses: 
                UC = so.minimize(objective, x0 = guess, method='SLSQP', 
                                 options={'maxiter': 2000, 'ftol': 1e-07,'disp': True}) 
                C = so.minimize(objective, x0 = guess, method='SLSQP',
                                options={'maxiter': 2000, 'ftol': 1e-07,'disp': True}, constraints=cons)
                            
                minimised_objective_at_guess_UC.append(UC)
                minimised_objective_at_guess_UC_fun.append(UC.fun) 
                minimised_objective_at_guess.append(C)
                minimised_objective_at_guess_fun.append(C.fun) 
           
            # define the desired solutions (those with the smallest minimised objective value) 
            sol_UC = minimised_objective_at_guess_UC[np.argmin(minimised_objective_at_guess_UC_fun)] 
            sol = minimised_objective_at_guess[np.argmin(minimised_objective_at_guess_fun)] 
            
            # > CONSTRAINED 
            if sol.success == True: 
                mT2dc_diff[0] = sol.fun - (alphaVal*mt2_W + (1-alphaVal)*mt2_t_11_22)
                mT2dc[0] = sol.fun
                mT2prime_W[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol.x)
                mT2prime_t[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol.x)
                sub_pT_sideA[0] = np.linalg.norm(sol.x)
                sub_pT_sideB[0] = np.linalg.norm(met[:2] - sol.x)
                sub_pT_min_over_met[0] = (min(np.linalg.norm(sol.x), np.linalg.norm(met[:2] - sol.x))/np.linalg.norm(met[:2]))
                alpha[0] = alphaVal 
                constraint_pT_cut[0] = pT_cut  
                calc_speed[0] = 0 
                success[0] = 1 
                
                tree.Fill() 
             
                if mT2_W < 1:
                    mT2prime_W_subC[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_alt.x) 
                    mT2prime_t_subC[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_alt.x)
                    constraint_pT_subcut[0] = pT_cut 
                    tree.Fill() 
            else: 
                sol_alt = so.minimize(objective, x0 = invis_sideA_array_guesses[np.argmin(minimised_objective_at_guess_fun)],
                                      method='COBYLA', options={'maxiter': 2000, 'ftol': 1e-07,'disp': True}, constraints=cons)
                
                mT2dc_diff[0] = sol_alt.fun - (alphaVal*mt2_W + (1-alphaVal)*mt2_t_11_22)
                mT2dc[0] = sol_alt.fun
                mT2prime_W[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_alt.x)
                mT2prime_t[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_alt.x)
                sub_pT_sideA[0] = np.linalg.norm(sol_alt.x)
                sub_pT_sideB[0] = np.linalg.norm(met[:2] - sol_alt.x)
                sub_pT_min_over_met[0] = (min(np.linalg.norm(sol_alt.x), np.linalg.norm(met[:2] - 
                                                                                        sol_alt.x))/np.linalg.norm(met[:2]))
                alpha[0] = alphaVal 
                constraint_pT_cut[0] = pT_cut  
                calc_speed[0] = 0 
                
                if mT2_W < 1:
                    mT2prime_W_subC[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_alt.x) 
                    mT2prime_t_subC[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_alt.x) 
                    constraint_pT_subcut[0] = pT_cut 
                    tree.Fill() 
                
                if sol_alt.success == True:
                    success[0] = 1
                else:
                    success[0] = 0 
                
                tree.Fill() 
                
            # > UNCONSTRAINED 
            mT2dc_diff[0] = sol_UC.fun - (alphaVal*mt2_W + (1-alphaVal)*mt2_t_11_22)
            mT2dc[0] = sol_UC.fun
            mT2prime_W[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x)
            mT2prime_t[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x)
            sub_pT_sideA[0] = np.linalg.norm(sol_UC.x)
            sub_pT_sideB[0] = np.linalg.norm(met[:2] - sol_UC.x)
            sub_pT_min_over_met[0] = (min(np.linalg.norm(sol_UC.x), np.linalg.norm(met[:2] - 
                                                                                        sol_UC.x))/np.linalg.norm(met[:2]))
            mT2_W[0] = mt2_W 
            mT2_t[0] = mt2_t 
            alpha[0] = alphaVal 
            constraint_pT_cut[0] = 0 
            calc_speed[0] = 0
            success[0] = 1 
            
            tree.Fill() 
 
            if mT2_W < 1:
                mT2prime_W_subC[0] = DC.get_alpha_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x) 
                mT2prime_t_subC[0] = DC.get_beta_term(vis_sideA_array, vis_sideB_array, met, sol_UC.x) 
                constraint_pT_subcut[0] = 0
                tree.Fill() 
                
                
                
##############################################
# Draw all histograms.
##############################################
c = ROOT.TCanvas()

#### INPUT TREE FILE #### 
h_ell1_pt.Draw("E") 
c.SaveAs(outDir+"h_ell1_PT.pdf")
h_ell1_E.Draw("E")
c.SaveAs(outDir+"h_ell1_E.pdf")

h_ell2_pt.Draw("E")
c.SaveAs(outDir+"h_ell2_PT.pdf")
h_ell2_E.Draw("E")
c.SaveAs(outDir+"h_ell2_E.pdf")

h_bjet1_E.Draw("E")
c.SaveAs(outDir+"h_bjet1_E.pdf")
h_bjet2_E.Draw("E")
c.SaveAs(outDir+"h_bjet2_E.pdf")

h_mT2_W.Draw("E")
c.SaveAs(outDir+"h_mT2_W.pdf")
h_mT2_t_11_22.Draw("E")
c.SaveAs(outDir+"h_mT2_t_11_22.pdf")
h_mT2_t_12_21.Draw("E")
c.SaveAs(outDir+"h_mT2_t_12_21.pdf")
h_mT2_t_min.Draw("E")
c.SaveAs(outDir+"h_mT2_t_min.pdf")

h_EtMiss.Draw("E") 
c.SaveAs(outDir+"h_EtMiss.pdf")
h_EtMiss_phi.Draw("E") 
c.SaveAs(outDir+"h_EtMiss_phi.pdf")

##############################################
# Save histograms to ROOT Output Files. 
##############################################
#### INPUT TREE FILE ####
h_ell1_pt.Write()
h_ell1_E.Write()

h_ell2_pt.Write()
h_ell2_E.Write()

h_bjet1_E.Write()
h_bjet2_E.Write()

h_mT2_W.Write()
h_mT2_t_11_22.Write()
h_mT2_t_12_21.Write()
h_mT2_t_min.Write()

h_EtMiss.Write()
h_EtMiss_phi.Write() 

f_outputRoot.Write() 
f_outputRoot.Close()