###############################################
# Mt2dc fitting analysis
##############################################

import ROOT
import numpy as np 

##############################################
# Define the input and output root files
##############################################
f_inputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_07_July_7_mt2dc_makeFinalPlots.root", "read")

outDir = "/Users/juliakim/Documents/styledPlotsOutputs/fit_functions/" 
parameterFile = open("/Users/juliakim/Documents/styledPlotsOutputs/fit_functions_parameters.txt",'w')
parameterFile.write('parameter 1, 2, 3, 4, 5, reduced chi-squared, derivative at 1 GeV, ∆x symm. drop-off \n')

##############################################
# Define constants and optimisation mode 
##############################################
m_W = 80 # GeV 
m_t = 173 # GeV

##############################################
# Get histograms for fitting 
##############################################
ROOT.gStyle.SetOptStat(1111)
c = ROOT.TCanvas()

#### Get histograms
h2_muT2dc = f_inputRoot.Get("h2_muT2dc").Clone()
h2_muT2dc_UC = f_inputRoot.Get("h2_muT2dc_UC").Clone()
h2_muT2prime_W = f_inputRoot.Get("h2_muT2prime_W").Clone()
h2_muT2prime_W_UC = f_inputRoot.Get("h2_muT2prime_W_UC").Clone()
h2_muT2prime_t = f_inputRoot.Get("h2_muT2prime_t").Clone()
h2_muT2prime_t_UC = f_inputRoot.Get("h2_muT2prime_t_UC").Clone()

#### Get number of alpha bins 
h_proj_x = h2_muT2dc.ProjectionX("h_proj_x", 0, -1)
num_alpha_bins = h_proj_x.GetNbinsX()

#### Loop over every alpha bin to create num_alpha_bins TH1 histograms per TH2 histogram
for i in range(num_alpha_bins):
    h_muT2dc_i = h2_muT2dc.ProjectionY("h_muT2dc_i", i, i)
    h_muT2dc_UC_i = h2_muT2dc_UC.ProjectionY("h_muT2dc_UC_i", i, i)
    h_muT2prime_W_i = h2_muT2prime_W.ProjectionY("h_muT2prime_W_i", i, i)
    h_muT2prime_W_UC_i = h2_muT2prime_W_UC.ProjectionY("h_muT2prime_W_UC_i", i, i)    
    h_muT2prime_t_i = h2_muT2prime_t.ProjectionY("h_muT2prime_t_i", i, i)    
    h_muT2prime_t_UC_i = h2_muT2prime_t_UC.ProjectionY("h_muT2prime_t_UC_i", i, i) 
    
    ### Fit every TH1 histogram & write data 
    f1 = ROOT.TF1("f1","([2]+[3]*x)*atan((x-[0])/[1]) + [4]", 0, 140) 
    f1.SetParameters(81., 13., -971., 2., 922.) 
    
    f2 = ROOT.TF1("f2", "([2] + [3]*x)*(x-[0])/sqrt([1]+(x-[0])**2) + [4]*x + 1", 0, 140)
    f2.SetParameters(75.,  434., -624.,  -13., 17.) 
    
    ## muT2dc
    # f1 
    h_muT2dc_i.Fit("f1", "0", "", 50, 140)
    f1.Draw("C")
    h_muT2dc_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2dc_" + str(i) + "_f1.pdf") 
    h_muT2dc_i_f1_data = np.around(np.array([f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2), f1.GetParameter(3),
                                             f1.GetParameter(4), f1.GetChisquare()/f1.GetNDF(), f1.Derivative(1), 
                                             (f1.GetX(f1(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2dc_' + str(i) + '_f1: \t' + str(h_muT2dc_i_f1_data) + '\n')
    
    # f2 
    h_muT2dc_i.Fit("f2", "0", "", 50, 140)
    f2.Draw("C")
    h_muT2dc_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2dc_" + str(i) + "_f2.pdf") 
    h_muT2dc_i_f2_data = np.around(np.array([f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2), f2.GetParameter(3),
                                             f2.GetParameter(4), f2.GetChisquare()/f2.GetNDF(), f2.Derivative(1), 
                                             (f2.GetX(f2(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2dc_' + str(i) + '_f2: \t' + str(h_muT2dc_i_f2_data) + '\n')
    
    ## muT2dc_UC 
    # f1 
    h_muT2dc_UC_i.Fit("f1", "0", "", 50, 140)
    f1.Draw("C")
    h_muT2dc_UC_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2dc_UC_" + str(i) + "_f1.pdf") 
    h_muT2dc_UC_i_f1_data = np.around(np.array([f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2), f1.GetParameter(3),
                                             f1.GetParameter(4), f1.GetChisquare()/f1.GetNDF(), f1.Derivative(1), 
                                             (f1.GetX(f1(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2dc_UC_' + str(i) + '_f1: \t' + str(h_muT2dc_UC_i_f1_data) + '\n')
    
    # f2 
    h_muT2dc_UC_i.Fit("f2", "0", "", 50, 140)
    f2.Draw("C")
    h_muT2dc_UC_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2dc_UC_" + str(i) + "_f2.pdf") 
    h_muT2dc_UC_i_f2_data = np.around(np.array([f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2), f2.GetParameter(3),
                                             f2.GetParameter(4), f2.GetChisquare()/f2.GetNDF(), f2.Derivative(1), 
                                             (f2.GetX(f2(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2dc_UC_' + str(i) + '_f2: \t' + str(h_muT2dc_UC_i_f2_data) + '\n')
    
    ## h_muT2prime_W
    # f1 
    h_muT2prime_W_i.Fit("f1", "0", "", 50, 140)
    f1.Draw("C")
    h_muT2prime_W_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2prime_W_" + str(i) + "_f1.pdf") 
    h_muT2prime_W_i_f1_data = np.around(np.array([f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2), f1.GetParameter(3),
                                             f1.GetParameter(4), f1.GetChisquare()/f1.GetNDF(), f1.Derivative(1), 
                                             (f1.GetX(f1(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2prime_W_' + str(i) + '_f1: \t' + str(h_muT2prime_W_i_f1_data) + '\n')
    
    # f2 
    h_muT2prime_W_i.Fit("f2", "0", "", 50, 140)
    f2.Draw("C")
    h_muT2prime_W_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2prime_W_" + str(i) + "_f2.pdf") 
    h_muT2prime_W_i_f1_data = np.around(np.array([f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2), f2.GetParameter(3),
                                             f2.GetParameter(4), f2.GetChisquare()/f2.GetNDF(), f2.Derivative(1), 
                                             (f2.GetX(f2(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2prime_W_' + str(i) + '_f2: \t' + str(h_muT2prime_W_i_f2_data) + '\n')
    
    ## h_muT2prime_W_UC
    # f1 
    h_muT2prime_W_UC_i.Fit("f1", "0", "", 50, 140)
    f1.Draw("C")
    h_muT2prime_W_UC_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2prime_W_UC_" + str(i) + "_f1.pdf") 
    h_muT2prime_W_UC_i_f1_data = np.around(np.array([f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2), 
                                                     f1.GetParameter(3), f1.GetParameter(4), f1.GetChisquare()/f1.GetNDF()
                                                     f1.Derivative(1), (f1.GetX(f1(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2prime_W_UC_' + str(i) + '_f1: \t' + str(h_muT2prime_W_UC_i_f1_data) + '\n')
    
    # f2 
    h_muT2prime_W_UC_i.Fit("f2", "0", "", 50, 140)
    f2.Draw("C")
    h_muT2prime_W_UC_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2prime_W_UC_" + str(i) + "_f2.pdf") 
    h_muT2prime_W_UC_i_f2_data = np.around(np.array([f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2),
                                                     f2.GetParameter(3),
                                                     f2.GetParameter(4), f2.GetChisquare()/f2.GetNDF(), f2.Derivative(1), 
                                                     (f2.GetX(f2(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2prime_W_UC_' + str(i) + '_f2: \t' + str(h_muT2prime_W_UC_i_f2_data ) + '\n')

    ## h_muT2prime_t
    # f1 
    h_muT2prime_t_i.Fit("f1", "0", "", 50, 140)
    f1.Draw("C")
    h_muT2prime_t_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2prime_t_" + str(i) + "_f1.pdf") 
    h_muT2prime_t_i_f1_data = np.around(np.array([f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2), f1.GetParameter(3),
                                             f1.GetParameter(4), f1.GetChisquare()/f1.GetNDF(), f1.Derivative(1), 
                                             (f1.GetX(f1(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2prime_t_' + str(i) + '_f1: \t' + str(h_muT2prime_t_i_f1_data) + '\n')
    
    # f2 
    h_muT2prime_t_i.Fit("f2", "0", "", 50, 140)
    f2.Draw("C")
    h_muT2prime_t_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2prime_t_" + str(i) + "_f2.pdf") 
    h_muT2prime_t_i_f1_data = np.around(np.array([f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2), f2.GetParameter(3),
                                             f2.GetParameter(4), f2.GetChisquare()/f2.GetNDF(), f2.Derivative(1), 
                                             (f2.GetX(f2(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2prime_t_' + str(i) + '_f2: \t' + str(h_muT2prime_t_i_f2_data) + '\n')
    

    ## h_muT2prime_t_UC
    # f1 
    h_muT2prime_t_UC_i.Fit("f1", "0", "", 50, 140)
    f1.Draw("C")
    h_muT2prime_t_UC_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2prime_t_UC_" + str(i) + "_f1.pdf") 
    h_muT2prime_t_UC_i_f1_data = np.around(np.array([f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2), 
                                                     f1.GetParameter(3), f1.GetParameter(4), f1.GetChisquare()/f1.GetNDF()
                                                     f1.Derivative(1), (f1.GetX(f1(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2prime_t_UC_' + str(i) + '_f1: \t' + str(h_muT2prime_t_UC_i_f1_data) + '\n')
    
    # f2 
    h_muT2prime_t_UC_i.Fit("f2", "0", "", 50, 140)
    f2.Draw("C")
    h_muT2prime_t_UC_i.Draw("E same") 
    c.SaveAs(outDir + "h_muT2prime_t_UC_" + str(i) + "_f2.pdf") 
    h_muT2prime_t_UC_i_f2_data = np.around(np.array([f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2),
                                                     f2.GetParameter(3), f2.GetParameter(4), f2.GetChisquare()/f2.GetNDF(), 
                                                     f2.Derivative(1), (f2.GetX(f2(50)/2, 50, 100) - 1)*2]), 3) 
    parameterFile.write('h_muT2prime_t_UC_' + str(i) + '_f2: \t' + str(h_muT2prime_t_UC_i_f2_data) + '\n')

parameterFile.close() 

   