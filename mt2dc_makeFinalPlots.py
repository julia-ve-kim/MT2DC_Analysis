#############################################################
# mt2dc final plotting make script
# Authors:
#    Ewan Hill
# 
# This script reads in histograms/plots and makes the final
# versions of them, some of which are overlays.
# 
# To run you need to make the styledPlotsOutputs/ directory
#############################################################

import ROOT

ROOT.gROOT.ForceStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()

##############################################
# Define the input and output root files
##############################################
inFileName = "/Users/juliakim/Documents/2022_06_June_15_truth_analysis_v01.root" 
inputFile = ROOT.TFile(inFileName, "read")
inFileName2 = '/Users/juliakim/Documents/2022_05_May_10_mt2dc_analysis_FOROVERLAY.root'
inputFile2 = ROOT.TFile(inFileName2, "read") 
outDir = "/Users/juliakim/Documents/styledPlotsOutputs/"  

##############################################
# Define constants 
##############################################
m_W = 80 # GeV 
m_t = 173 # GeV

##############################################
# Mt2dc Overlay Code 
##############################################
cbf_Black         = ROOT.TColor.GetFreeColorIndex()
colour1           = ROOT.TColor(cbf_Black,            0.00, 0.00, 0.00, "cbf_Black")
cbf_Orange        = ROOT.TColor.GetFreeColorIndex()
colour2           =  ROOT.TColor(cbf_Orange,          0.90, 0.60, 0.00, "cbf_Orange")
cbf_SkyBlue       = ROOT.TColor.GetFreeColorIndex()
colour3           =  ROOT.TColor(cbf_SkyBlue,         0.35, 0.70, 0.90, "cbf_SkyBlue")
cbf_BluishGreen   = ROOT.TColor.GetFreeColorIndex()
colour4           =  ROOT.TColor(cbf_BluishGreen,     0.00, 0.60, 0.50, "cbf_BluishGreen")
cbf_Yellow        = ROOT.TColor.GetFreeColorIndex()
colour5           =  ROOT.TColor(cbf_Yellow,          0.95, 0.90, 0.25, "cbf_Yellow")
cbf_Blue          = ROOT.TColor.GetFreeColorIndex()
colour6           =  ROOT.TColor(cbf_Blue,            0.00, 0.45, 0.70, "cbf_Blue")
cbf_Vermilion     = ROOT.TColor.GetFreeColorIndex()
colour7           =  ROOT.TColor(cbf_Vermilion,       0.80, 0.40, 0.00, "cbf_Vermilion")
cbf_ReddishPurple = ROOT.TColor.GetFreeColorIndex()
colour8           =  ROOT.TColor(cbf_ReddishPurple,   0.80, 0.60, 0.70, "cbf_ReddishPurple")

# c1 = new TCanvas("c1", "", 0, 0, 1000, 700)
c1 = ROOT.TCanvas("c1", "")
ROOT.gROOT.ForceStyle() 
#c1.SetLogy() # use 60 to set logarithmic y-axis scale 

plotCounter = 0
plotCounter += 1
print("plotCounter = ", plotCounter)

histLabel__mtW_mt2 = "#Sigma = max(#mbox{mT_{1}^{true}}, #mbox{mT_{2}^{true}})" 
histLabel__mtW_mt2prime = "#mbox{m_{T2}}(W)" 
inputHist__mtW_mt2 = "h_max_mT_true"  
inputHist__mtW_mt2prime = "h_mT2_W" # h_mT2dc_alpha_1 (equivalent to) 
xAxisDescription__mtW = "Transverse mass estimator [GeV]" 
yAxisDescription__mtW = "Number of events / 1 GeV" 
outFile__mtW = "mT2_W_max_mT_true_overlay" 
title__mtW = ""
ratioMin = 0.1
ratioMax = 20

h_v1 = inputFile.Get(inputHist__mtW_mt2).Clone()
h_v1.Scale(1./h_v1.Integral(), "width")
h_v1.SetTitle("")
h_v2 = inputFile2.Get(inputHist__mtW_mt2prime).Clone()
h_v2.Scale(1./h_v2.Integral(), "width")
h_v2.SetTitle("")

h_v1.GetXaxis().SetTitle(xAxisDescription__mtW)
h_v1.GetYaxis().SetTitle(yAxisDescription__mtW)
h_v2.GetXaxis().SetTitle(xAxisDescription__mtW)
h_v2.GetYaxis().SetTitle(yAxisDescription__mtW)

numEntries_v1 = h_v1.Integral()
numEntries_v2 = h_v2.Integral()
print("numEntries_v1 = ", numEntries_v1, ", numEntries_v2 = ", numEntries_v2)
maxBinYvalue_h_v1 = h_v1.GetMaximum()
maxBinYvalue_h_v2 = h_v2.GetMaximum()
underflow_v1 = h_v1.GetBinContent( 0 )
underflow_v2 = h_v2.GetBinContent( 0 )
overflow_v1 = h_v1.GetBinContent( h_v1.GetXaxis().GetNbins() + 1 )
overflow_v2 = h_v2.GetBinContent( h_v2.GetXaxis().GetNbins() + 1 )

maxBinValue = 1.25 * max(maxBinYvalue_h_v1, maxBinYvalue_h_v2)   

# Rebinning
h_v1.Rebin(1)
h_v2.Rebin(1)

# Use this as part of shifting the y axis title position and not having it overlap the rest of the plot
c1.SetLeftMargin(0.1)
c1.SetTicks(1, 1)   # ticks on top part of ratio plot

# SetMoreLogLabels()   # how to make this work in pyroot ?

# Use this for shifting the y axis title position
h_v1.GetYaxis().SetTitleOffset(1.5)

# Make a ratio plot
rp = ROOT.TRatioPlot(h_v1, h_v2)   # ratio of histograms in the bottom panel
# dp = ROOT.TRatioPlot(h_v2, h_v1, "diff")   # difference of histograms in the bottom pannel
c1.SetTicks(1, 1)   # ticks on top part of ratio plot
rp.Draw()

# Format the ratio plot in the bottom panel
rp.GetLowerRefGraph().SetMinimum(ratioMin)
rp.GetLowerRefGraph().SetMaximum(ratioMax)
rp.GetLowerRefXaxis().SetTitle(xAxisDescription__mtW)
rp.GetLowerRefYaxis().SetTitle("#Sigma / #mbox{m_{T2}}(W)") 

# Format the overlayed histograms in the top panel
rp.GetUpperPad().cd()

h_v1.SetLineColor(cbf_Black)
h_v1.SetFillColor(cbf_SkyBlue)

h_v2.SetMarkerColor(cbf_Black)
h_v2.SetLineColor(cbf_Black)
h_v2.SetMarkerStyle(20)
h_v2.SetMarkerSize(1.0)

rp.GetUpperPad().Update()

# Plot title - implemented using TText
t = ROOT.TText(.18,.92,title__mtW)
t.SetTextFont(43)
t.SetTextSize(12)
t.SetNDC()
t.Draw()

# Create a legend to define the different histograms for the reader
leg = ROOT.TLegend(.60, .60, .75, .80)
leg.SetFillColor(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.038)

leg.AddEntry(h_v1, histLabel__mtW_mt2, "f")
leg.AddEntry(h_v2, histLabel__mtW_mt2prime, "PE1")

leg.Draw()

c1.Update()   # Final update, just in case.

# Output the final plot in various forms.  Want several formats
#    since we will want pdf as the final version for presenting
#    but may need to make modifications, which can sometimes
#    be most easily done manually.
c1.SaveAs(outDir + outFile__mtW + ".pdf")
c1.SaveAs(outDir + outFile__mtW + ".C")
c1.SaveAs(outDir + outFile__mtW + ".root")

inputFile.Close()

print("Finished")