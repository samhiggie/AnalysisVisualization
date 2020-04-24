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
    from utils.Processes import HAA_processes
    name = h.GetName()
    if "a15" in name or "a20" in name or "a25" in name or "a30" in name or "a35" in name or "a40" in name: 
        h.SetLineColor(ROOT.kRed)
        h.SetFillColor(0)
        h.SetLineWidth(4)
    if "ZZ" in name: 
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["ZZ"].color[0]))
    if "ZTT" in name: 
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["DYJetsToLL"].color[0]))
    if "TT" in name: 
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["TTJets_1LTbar"].color[0]))
    if "Z" in name: 
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["DYJetsToLL"].color[0]))
    if "W" in name: 
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["WJetsToLNu"].color[0]))
    if "EWK" in name: 
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["EWKWMinus2Jets"].color[0]))
    if "ST" in name: 
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["ST_s"].color[0]))
    if "TTL" in name: 
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["TTJets_1LTbar"].color[0]))
    if "ZTL" in name: 
        h.SetFillColor(ROOT.TColor.GetColor(HAA_processes["DYJetsToLL"].color[0]))
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
    from utils.Categories import allcats
    from utils.Processes import HAA_processes
    import argparse

    parser = argparse.ArgumentParser(description="make full plots from root files containing histograms")
    #parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
    parser.add_argument("-o",  "--output", default="",  help="postfix string")
    args = parser.parse_args()

    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)


    #style1=GetStyleHtt()
    #style1.cd()
    #c=MakeCanvas("asdf","asdf",750,750)


    #myfile=ROOT.TFile(args.input,"read")  
    myfile=ROOT.TFile("AMass.root","read")   
    
    #gathering all the categories as a list of files
    #for cat in allcats:
    catnames = HAA_Inc_mmmt.name
    newvars=[]
    variabledic={}
    for newvar in  HAA_Inc_mmmt.newvariables.keys():
        newvars.append([newvar])
        variabledic[newvar]=[newvar,HAA_Inc_mmmt.newvariables[newvar][1],HAA_Inc_mmmt.newvariables[newvar][3],HAA_Inc_mmmt.newvariables[newvar][4]]
    variables = [] 
    for varHandle in HAA_Inc_mmmt.vars.keys():
        variables.append([varHandle])
        variabledic[varHandle]=HAA_Inc_mmmt.vars[varHandle]
    variables = variables + newvars
    #print variables
    filelist = {}
    processes = []
    
    for var in variables:
        filelist[var[0]] = var[0]+".root"

    for pro in HAA_processes.values():
        for subpro in pro.cuts.keys():
            processes.append(subpro)
    processes=list(dict.fromkeys(processes))
    #print "Will plot the following physics processes ",processes


    for ivar,var in enumerate(filelist.keys()):
    #for ivar,var in enumerate(cat.vars.keys()):
        print "Working on file ",filelist[var],"  for variable   ",var
        file = ROOT.TFile("out/"+filelist[var],"read")
        #bins=HAA_Inc_mmmt.binning + HAA
         
        #for catname in enumerate():
        cati = 0
        for categoryKey in file.GetListOfKeys():
            category = categoryKey.ReadObj()
            print "Working on category ",category.GetName()
            try:
                os.mkdir("outplots_"+args.output+"_"+str(category.GetName()))
            except:
                print "dir prob exists"
            
            if var=="AMass":
                fileout = open("outplots_"+args.output+"_"+str(category.GetName())+"/"+str(allcats[cati].name)+"_info.txt","w")  
                fileout.write("Working on category "+allcats[cati].name+"\n")

            #signal 
            hSignal = category.Get("a40")
            
            #reducible
            #hBackground = ROOT.THStack("reducible background","")
            hBackground = category.Get("Z")
            hBackground.Add(category.Get("W"))
            hBackground.Add(category.Get("TT"))
            hBackground.Add(category.Get("ST"))
            hBackground.Add(category.Get("EWK"))
            
            #irreducible
            #hirBackground = ROOT.THStack("irreducible background","")
            #hirBackground = ROOT.TH1D()
            hirBackground = category.Get("ZZ")
            hirBackground.Add(category.Get("ZHToTauTau"))

            #triple coupling 
            h3alphBackground = category.Get("ttZ")
            h3alphBackground.Add(category.Get("WWZ"))
            h3alphBackground.Add(category.Get("WZZ"))
            h3alphBackground.Add(category.Get("ZZZ"))
            h3alphBackground.Add(category.Get("HZJ"))
            h3alphBackground.Add(category.Get("GluGluToContinToZZTo2mu2tau"))

            #rare 
            hrareBackground = category.Get("rare")
            #hrareBackground = category.Get("TTGamma_Hadr")
            #hrareBackground.Add(category.Get("TTGamma_1LT"))
            #hrareBackground.Add(category.Get("TTGamma_1LTbar"))
            #hrareBackground.Add(category.Get("TTGamma_2L"))
            #hrareBackground.Add(category.Get("TTHH"))
            #hrareBackground.Add(category.Get("TTTJ"))
            #hrareBackground.Add(category.Get("TTTT"))
            #hrareBackground.Add(category.Get("TTTW"))
            #hrareBackground.Add(category.Get("TTW"))
            #hrareBackground.Add(category.Get("TTWH"))
            #hrareBackground.Add(category.Get("TTWJetsToQQ"))
            #hrareBackground.Add(category.Get("TTWW"))
            #hrareBackground.Add(category.Get("TTWZ"))
            #hrareBackground.Add(category.Get("TTZH"))
            #hrareBackground.Add(category.Get("TTZZ"))
            #hrareBackground.Add(category.Get("WW"))
            #hrareBackground.Add(category.Get("WWext"))
            #hrareBackground.Add(category.Get("WZext"))

            #data
            hData = category.Get("data_obs")


            
            #print "Working on category ... ",category.GetName()
            category.cd()
            histolist=[]

            hBackground.SetLineColor(1)
            hirBackground.SetLineColor(1)
            hrareBackground.SetLineColor(1)
            hBackground.SetFillStyle(1001)
            hirBackground.SetFillStyle(1001)
            hrareBackground.SetFillStyle(1001)
            hBackground.SetFillColor(ROOT.TColor.GetColor("#CF8AC8"))
            hirBackground.SetFillColor(ROOT.TColor.GetColor("#13E2FE"))
            #h3alphBackground.SetFillColor(ROOT.TColor.GetColor("#65E114"))
            hrareBackground.SetFillColor(ROOT.TColor.GetColor("#65E114"))
            #hrareBackground.SetFillColor(ROOT.TColor.GetColor("#FF6600"))
            h3alphBackground.SetFillColor(ROOT.TColor.GetColor("#FF6600"))
                        
            #hSignal.SetTitle()
            hBackground.SetTitle("reducible")
            hirBackground.SetTitle("irreducible")
            h3alphBackground.SetTitle("triple #alpha")
            hrareBackground.SetTitle("rare")

            #for ratio
            hbkgtot = hBackground.Clone()
            #hbkgtot = h3alphBackground.Clone()
            hbkgtot.Add(hirBackground)
            hbkgtot.Add(hrareBackground)

            hBkgTot = ROOT.THStack()
            hBkgTot.Add(hBackground)
            hBkgTot.Add(hirBackground)
            hBkgTot.Add(hrareBackground)
            #hBkgTot.Add(h3alphBackground)


            H = 600
            W = 600

            H_ref = 600
            W_ref = 600

           
            T = 0.08*H_ref
            B = 0.12*H_ref 
            L = 0.16*W_ref
            R = 0.04*W_ref

            
            B_ratio = 0.1*H_ref
            T_ratio = 0.03*H_ref

            B_ratio_label = 0.3*H_ref

            c=ROOT.TCanvas("canvas","",0,0,600,600)

            doRatio=True

            if not doRatio:
                c.SetLeftMargin     (L/W)
                c.SetRightMargin    (R/W)
                c.SetTopMargin      (T/H)
                c.SetBottomMargin   (B/H)
            c.cd()

            pad1 = ROOT.TPad("pad1","pad1",0.0016,0.291,1.0,1.0)
            #pad2 = ROOT.TPad("pad2","pad2",0.45,0.35,0.92,0.95)
            pad2 = ROOT.TPad("pad2","The lower pad",0,0,1.0,0.29)

            if(doRatio):

                pad1.SetTicks(0,0)
                pad1.SetLeftMargin(L/W)
                pad1.SetRightMargin(R/W)
                pad1.SetTopMargin(T/H)
                pad1.SetBottomMargin(B_ratio/H) 
                pad1.SetFillColor(0)
                pad1.SetBottomMargin(0)


                pad2.SetLeftMargin(L/W)
                pad2.SetRightMargin(R/W)
                pad2.SetTopMargin(T_ratio/H)
                pad2.SetTopMargin(0.007)
                pad2.SetBottomMargin(B_ratio_label/H)
                pad2.SetGridy(1)
                pad2.SetFillColor(4000)

            else:
                pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
                pad1.SetLeftMargin     (L/W)
                pad1.SetRightMargin    (R/W)
                pad1.SetTopMargin      (T/H)
                pad1.SetBottomMargin   (B/H)



            pad1.Draw()
            pad1.cd()
            lumi=add_lumi()
            cms=add_CMS()
            pre=add_Preliminary()
            xR=0.65
            l=ROOT.TLegend(xR,0.55,xR+0.28,0.9);
            l.AddEntry(hData)
            l.AddEntry(hSignal)
            l.AddEntry(hBackground)
            l.AddEntry(hirBackground)
            l.AddEntry(hrareBackground)
            #l.AddEntry(h3alphBackground)

            #print "signal entries   ",hSignal.GetEntries()
            hData.SetTitle("")
            hData.SetMarkerStyle(8)
            hData.Draw("ep")
            hBkgTot.Draw("HIST same")
            hSignal.Draw("same")
            hBkgTot.GetXaxis().SetTitle(variabledic[var][3]+variabledic[var][2])
            hBkgTot.GetYaxis().SetTitle("Events")
            hData.GetXaxis().SetTitle(variabledic[var][3]+variabledic[var][2])
            hData.GetYaxis().SetTitle("Events")
            hData.Draw("ep same")
            #hBackground.Draw("same")
            #hirBackground.Draw("same")
            #for hist in histolist:
            #    hist.Draw("same")
            lumi.Draw()
            cms.Draw()
            pre.Draw()
            l.Draw()            

            #TPad 2 for ratio
            c.cd()
            hData2=hData.Clone()
            pad2.Draw()
            pad2.cd()
            hData2.Divide(hbkgtot)
            hData2.SetMarkerStyle(20)
            hData2.SetTitleSize  (0.12,"Y")
            hData2.SetTitleOffset(0.40,"Y")
            hData2.SetTitleSize  (0.12,"X")
            hData2.SetLabelSize  (0.10,"X")
            hData2.SetLabelSize  (0.08,"Y")
            #hData2.GetYaxis().SetRangeUser(0.62,1.38)
            hData2.GetYaxis().SetRangeUser(0.3,2.2)
            hData2.GetYaxis().SetNdivisions(305)
            hData2.GetYaxis().SetTitle("Obs/Exp   ")
            hData2.Draw("P");
            pad2.Draw()

            #print "with cuts ",allcats[cati].cuts
            #print "data entries ",hData.GetEntries()
            #print "background entries ",hBackground.GetEntries()
            if var=="AMass":
                for key in allcats[cati].cuts.keys():
                    for cut in allcats[cati].cuts[key]:
                        fileout.write(str(cut))
                    fileout.write("\n")
                fileout.write("data entries "+str(hData.Integral())+"\n")
                fileout.write("background entries "+str(hBackground.Integral())+"\n")
            
            c.SaveAs("outplots_"+args.output+"_"+str(category.GetName())+"/"+var+"_"+str(category.GetName())+".png")
            cati = cati +1

            hBackground.Delete()
            hirBackground.Delete()
            h3alphBackground.Delete()
            hData.Delete()
            hData2.Delete()
            hBkgTot.Delete()


