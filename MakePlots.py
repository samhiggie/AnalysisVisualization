#!/usr/bin/env python
import ROOT
import re
import math
from array import array
#from collections import OrderedDict
#import varCfgPlotter
import argparse
import os
#from HttStyles import GetStyleHtt
#from HttStyles import MakeCanvas

def add_lumi():
    lowX=0.65
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.04)
    lumi.SetTextFont (   42 )
    lumi.AddText("35.9 fb^{-1} (13 TeV)")
    return lumi

def add_CMS():
    lowX=0.17
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.06)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS")
    return lumi

def add_Preliminary():
    lowX=0.45
    #lowY=0.690
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(52)
    lumi.SetTextSize(0.04)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("Preliminary")
    return lumi

def make_legend():
        output = ROOT.TLegend(0.165, 0.45, 0.350, 0.75, "", "brNDC")
        output.SetLineWidth(0)
        output.SetLineStyle(0)
        output.SetFillStyle(0)
        output.SetFillColor(0)
        output.SetBorderSize(0)
        output.SetTextFont(62)
        return output

def make_legend_inset():
        output = ROOT.TLegend(0.6, 0.5, 0.98, 0.85, "", "brNDC")
        output.SetLineWidth(0)
        output.SetLineStyle(0)
        output.SetFillStyle(0)
        output.SetFillColor(0)
        output.SetBorderSize(0)
        output.SetTextFont(62)
        return output

def add_categ(text):
       categ  = ROOT.TPaveText(0.6, 0.2+0.013, 0.89, 0.40+0.155, "NDC")
       categ.SetBorderSize(   0 )
       categ.SetFillStyle(    0 )
       categ.SetTextAlign(   12 )
       categ.SetTextSize ( 0.05 )
       categ.SetTextColor(    1 )
       categ.SetTextFont (   41 )
       categ.AddText(text)
       return categ

def normalize_hist(h):
    k=h.Clone()
    for j in range(1,h.GetSize()-1):
        k.SetBinContent(j,k.GetBinContent(j)/k.GetBinWidth(j))
        k.SetBinError(j,k.GetBinError(j)/k.GetBinWidth(j))
    return k

def assignColor(h): 
    name = h.GetName()
    if "a15" in name or "a20" in name or "a25" in name or "a30" in name or "a35" in name or "a40" in name: 
        h.SetFillColor(ROOT.kRed)
    if "ZL" in name: 
        h.SetFillColor()
    if "ZTT" in name: 
        h.SetFillColor()
    if "TTT" in name: 
        h.SetFillColor()
    if "TTL" in name: 
        h.SetFillColor()
    if "ZTL" in name: 
        h.SetFillColor()
    if "jetFakes" in name: 
        h.SetFillColor()
    return h 

def makeHisto(h,hS,hB,hD):
    #k=h.Clone()
    name = h.GetName()
    print "setting attributes for ",name
    if "a15" in name or "a20" in name or "a25" in name or "a30" in name or "a35" in name or "a40" in name: 
        h = assignColor(h)
        hS.Add(h)
    if "ZL" in name or "ZTT" in name or "TTT" in name or "TTL" in name or "ZTL" in name or "jetFakes" in name: 
        h = assignColor(h)
        hB.Add(h)
        
     


    return h,hS,hB,hD








if __name__ == "__main__":

    #Change these to input category files ?? via parser?
    from utils.Categories import HAA_Inc_mmmt
    from utils.Processes import HAA

    #parser = argsparse.ArgumentParser(description="make full plots from root files containing histograms")
    #parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
    #parser.add_argument("-i",  "--input", default="mvis.root",  help="input file containing histograms")
    #args = parser.parse_args()

    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)


    #style1=GetStyleHtt()
    #style1.cd()
    #c=MakeCanvas("asdf","asdf",750,750)


    #myfile=ROOT.TFile(args.input,"read")  
    myfile=ROOT.TFile("mmmt_inclusive_AMass.root","read")   
    
    catnames = HAA_Inc_mmmt.name
    variables = HAA_Inc_mmmt.variables  
    filelist = {}
    processes = []
    
    for var in variables:
        filelist[var[0]] = catnames+"_"+var[0]+".root"
    for pro in HAA.values():
        processes.append(pro[0])
    print "Will plot the following physics processes ",processes


    for var in filelist.keys():
        print "Working on file ",filelist[var],"  for variable   ",var
        file = ROOT.TFile(filelist[var],"read")
        for categoryKey in file.GetListOfKeys():
            category = categoryKey.ReadObj()
            if category:
                hSignal = ROOT.THStack("signal","")
                hBackground = ROOT.THStack("background","")
                hData = ROOT.TH1D()
                print "Working on category ... ",category.GetName()
                category.cd()
                histolist=[]
                for histokey in category.GetListOfKeys():
                    histoOriginal = histokey.ReadObj()
                    histo = histoOriginal.Clone()
                    print histo.GetName()
                    histolist.append(histo)
                    name = histo.GetName()
                    
                    if (name in processes):
                        #make classification here ?? or use the processes from Processes.py?
                        #histo,hSignal,hBackground,hData = makeHisto(histo,hSignal,hBackground,hData) 
                        if "a15" in name or "a20" in name or "a25" in name or "a30" in name or "a35" in name or "a40" in name: 
                            histo = assignColor(histo)
                            print histo.GetEntries()
                            hSignal.Add(histo)
                        if "ZL" in name or "ZTT" in name or "TTT" in name or "TTL" in name or "ZTL" in name or "jetFakes" in name: 
                            histo = assignColor(histo)
                            hBackground.Add(histo)
                        
                    #else (histo.IsA().InheritsFrom("TH1")):
                    #    cother=ROOT.TCanvas("canvasOther","",0,0,600,600)
                    #    cother.cd()
                    #    histo.Draw()
                    #    cother.SaveAs("other_"+histo.GetName()+".png")
                    #    cother.Delete()
                    #histo.Delete()
                c=ROOT.TCanvas("canvas","",0,0,600,600)
                pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
                pad2 = ROOT.TPad("pad2","pad2",0.45,0.35,0.92,0.95)
                pad2 = ROOT.TPad("pad2","The lower pad",0,0,1.0,0.29)
                c.cd()
                pad1.Draw()
                pad1.cd()
                lumi=add_lumi()
                cms=add_CMS()
                pre=add_Preliminary()
                #print "signal entries   ",hSignal.GetEntries()
                hSignal.Draw()
                hSignal.GetXaxis().SetTitle(str(var))
                #for hist in histolist:
                #    hist.Draw("same")
                lumi.Draw()
                cms.Draw()
                pre.Draw()
        
                c.SaveAs(var+".png")
                c.Delete()
                pad1.Delete()
                pad2.Delete()
                #hBackground.Draw()
                #hData.Draw()
                #pad2.cd()



