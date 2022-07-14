###############################################
# Mt2dc fitting analysis
##############################################

import ROOT
import numpy as np 

##############################################
# Define the input and output root files
##############################################
inputFile = ROOT.TFile.Open("/Users/juliakim/Documents/2022_05_May_10_mt2dc_analysis_v01.root", "read")

outDir = "/Users/juliakim/Documents/styledPlotsOutputs/" 
parameterFile = open("/Users/juliakim/Documents/styledPlotsOutputs/fit_functions_parameters.txt",'w')

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
    
inputHist = "h_muT2dc_alpha_0"  # description: alpha = 1, constrained = lower bound of 20, 50 to 100 
h = inputFile.Get(inputHist).Clone()

# FIT 0 
f0 = ROOT.TF1("f0","([2]+[3]*x)*atan((x-[0])/[1]) + [4]", 0, 140) 
f0.SetParameters(81., 13., -971., 2., 922.) # (80, 10, -800, 0, 800) 
h.Fit("f0", "0", "", 50, 140)
f0.Draw("C")
h.Draw("E same") 
c.SaveAs(outDir + "f0.pdf") 

# FIT 1 
f1 = ROOT.TF1("f1","([2]+[3]*x)*atan((x-[0])/[1]) + [4]*x + 1", 0, 140) 
f1.SetParameters(76., 15., -731., -6., 15.) # (80, 10, -800, -10, 10) 
h.Fit("f1", "0", "", 50, 140)
f1.Draw("C")
h.Draw("E same") 
c.SaveAs(outDir + "f1.pdf") 

# FIT 2
f2 = ROOT.TF1("f2", "([2] + [3]*x)*(x-[0])/sqrt([1]+(x-[0])**2) + [4]", 0, 140)
f2.SetParameters(80, 321, -1200,  2,  958) # (80, 300, -800, 0, 800)
h.Fit("f2", "0", "", 50, 140)
f2.Draw("C")
h.Draw("E same")
c.SaveAs(outDir + "f2.pdf")

# FIT 3 
f3 = ROOT.TF1("f3", "([2] + [3]*x)*(x-[0])/sqrt([1]+(x-[0])**2) + [4]*x + 1", 0, 140)
f3.SetParameters(75.,  434., -624.,  -13., 17.) # (80, 300, -800, -10, 10) 
h.Fit("f3", "0", "", 50, 140)
f3.Draw("C")
h.Draw("E same")
c.SaveAs(outDir + "f3.pdf")

# FIT 4 
f4 = ROOT.TF1("f4","[0]*atan([1]*(x-[2]))+[3]", 0, 140) 
f4.SetParameters(-674., 0., 80.,  958.) # (-800, 0.1, 80, 800)
h.Fit("f4", "0", "", 50, 140)
f4.Draw("C")
h.Draw("E same") 
c.SaveAs(outDir + "f4.pdf") 

##############################################
# Drop-off analysis 
##############################################
# WEIGHTED AVERAGE
# see main code, 'mt2dc_analysis_v04.' 
    
# DERIVATIVE 
f0_slope = f0.Derivative(80) 
f1_slope = f1.Derivative(80)
f2_slope = f2.Derivative(80)
f3_slope = f3.Derivative(80)
f4_slope = f4.Derivative(80) 

# SYMMETRIC STEP-OUT  
f0_initial_y = f0(50)  
f0_drop_off_x = f0.GetX(f0_initial_y/2, 50, 100) # x value by which f1 drops off by half 
f0_delta_x = (f0_drop_off_x - 80)*2 # symmetric interval about 80 GeV 

f1_initial_y = f1(50)  
f1_drop_off_x = f1.GetX(f1_initial_y/2, 50, 100) # x value by which f1 drops off by half 
f1_delta_x = (f1_drop_off_x - 80)*2 # symmetric interval about 80 GeV 

f2_initial_y = f2(50)  
f2_drop_off_x = f2.GetX(f2_initial_y/2, 50, 100) 
f2_delta_x = (f2_drop_off_x - 80)*2

f3_initial_y = f3(50)  
f3_drop_off_x = f3.GetX(f3_initial_y/2, 50, 100) 
f3_delta_x = (f3_drop_off_x - 80)*2

f4_initial_y = f4(50)  
f4_drop_off_x = f4.GetX(f4_initial_y/2, 50, 100) 
f4_delta_x = (f4_drop_off_x - 80)*2

##############################################
# Save all data to output files 
##############################################

f0_parameters = np.around(np.array([f0.GetParameter(0), f0.GetParameter(1), f0.GetParameter(2), f0.GetParameter(3), f0.GetParameter(4), f0.GetChisquare()/f0.GetNDF(), f0_slope, f0_delta_x]), 1) 
f1_parameters = np.around(np.array([f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2), f1.GetParameter(3), f1.GetParameter(4), f1.GetChisquare()/f1.GetNDF(), f1_slope, f1_delta_x]), 1) 
f2_parameters = np.around(np.array([f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2), f2.GetParameter(3), f2.GetParameter(4), f2.GetChisquare()/f2.GetNDF(), f2_slope, f2_delta_x]), 1) 
f3_parameters = np.around(np.array([f3.GetParameter(0), f3.GetParameter(1), f3.GetParameter(2), f3.GetParameter(3), f3.GetParameter(4), f3.GetChisquare()/f3.GetNDF(), f3_slope, f3_delta_x]), 1) 
f4_parameters = np.around(np.array([f4.GetParameter(0), f4.GetParameter(1), f4.GetParameter(2), f4.GetParameter(3), f4.GetChisquare()/f4.GetNDF(), f4_slope, f4_delta_x]), 1)

parameterFile.write('h_mT2dc_alpha_1 [constraint = 20 GeV] \n') 
parameterFile.write('parameter 1, 2, 3, 4, 5, reduced chi-squared, derivative at 80 GeV, ∆x symm. drop-off \n')
parameterFile.write('f0_parameters: \t' + str(f0_parameters) + '\n')
parameterFile.write('f1_parameters: \t' + str(f1_parameters) + '\n')
parameterFile.write('f2_parameters: \t' + str(f2_parameters) + '\n')
parameterFile.write('f3_parameters: \t' + str(f3_parameters) + '\n')
parameterFile.write('f4_parameters: \t' + str(f4_parameters) + '\n')

parameterFile.close() 