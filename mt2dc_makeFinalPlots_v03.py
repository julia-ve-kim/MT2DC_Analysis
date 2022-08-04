###############################################
# Mt2dc fitting analysis
##############################################

import ROOT
import numpy as np 

##############################################
# Define the input and output root files
##############################################
alphaList = np.linspace(0, 1, 20) 
 
mT2prime_W_UC_f1_Parameter0 = np.loadtxt("/Users/juliakim/Documents/styledPlotsOutputs/mT2prime_W_UC_f1_Parameter0.txt", unpack=True, skiprows = 1, delimiter=',', usecols = 0)
mT2prime_W_UC_f1_Parameter0_ucty = np.loadtxt("/Users/juliakim/Documents/styledPlotsOutputs/mT2prime_W_UC_f1_Parameter0.txt", unpack=True, skiprows = 1, delimiter=',', usecols = 1)
mT2prime_W_UC_f2_Parameter0 = np.loadtxt("/Users/juliakim/Documents/styledPlotsOutputs/mT2prime_W_UC_f2_Parameter0.txt", unpack=True, skiprows = 1, delimiter=',', usecols=0)
mT2prime_W_UC_f2_Parameter0_ucty = np.loadtxt("/Users/juliakim/Documents/styledPlotsOutputs/mT2prime_W_UC_f2_Parameter0.txt", unpack=True, skiprows = 1, delimiter=',', usecols=1)

outDir = "/Users/juliakim/Documents/styledPlotsOutputs/"

c = ROOT.TCanvas()
gr = ROOT.TGraph(20, alphaList, mT2prime_W_UC_f1_Parameter0) 
gr.Draw("AC*") 
gr.SetTitle("mass_fit1(W) vs. #alpha; #alpha; mass_fit1(W) [GeV]") 
c.SaveAs(outDir + "mT2prime_W_UC_f1_Parameter0_vs_alpha.pdf")

gr2 = ROOT.TGraph(20, alphaList, mT2prime_W_UC_f1_Parameter0_ucty) 
gr2.Draw("AC*") 
gr2.SetTitle("mass_fit1_ucty(W) vs. #alpha; #alpha; mass_fit1_ucty(W) [GeV]") 
c.SaveAs(outDir + "mT2prime_W_UC_f1_Parameter0_ucty_vs_alpha.pdf") 
    
gr3 = ROOT.TGraph(20, alphaList, mT2prime_W_UC_f2_Parameter0)
gr3.Draw("AC*")
gr3.SetTitle("mass_fit2(W) vs. #alpha; #alpha; mass_fit2(W) [GeV]") 
c.SaveAs(outDir + "mT2prime_W_UC_f2_Parameter0_vs_alpha.pdf") 

gr4 = ROOT.TGraph(20, alphaList, mT2prime_W_UC_f2_Parameter0_ucty)
gr4.Draw("AC*")
gr4.SetTitle("mass_fit2_ucty(W) vs. #alpha; #alpha; mass_fit2_ucty(W) [GeV]") 
c.SaveAs(outDir + "mT2prime_W_UC_f2_Parameter0_ucty_vs_alpha.pdf") 
  
gr5 = ROOT.TGraph(20, alphaList, np.abs(mT2prime_W_UC_f1_Parameter0-mT2prime_W_UC_f2_Parameter0))
gr5.Draw("AC*")
gr5.SetTitle("|mass_fit1(W) - mass_fit2(W)| vs. #alpha; #alpha; |mass_fit1(W) - mass_fit2(W)| [GeV]") 
c.SaveAs(outDir + "mT2prime_W_UC_f1_f2_vs_alpha.pdf") 
