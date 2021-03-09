#!/usr/bin/env python
import ROOT
import re
import math
from array import array
import numpy as np
#from collections import OrderedDict
#import varCfgPlotter
import argparse
import os
import io
import yaml
#from HttStyles import GetStyleHtt
#from HttStyles import MakeCanvas
import uproot
import root_numpy

def add_lumi(year):
    lowX=0.65
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.04)
    lumi.SetTextFont (   42 )
    lumi.AddText(str(year)+" 35.9 fb^{-1} (13 TeV)")
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

def add_Preliminary(channel="mmmt"):
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
    lumi.AddText("Preliminary "+channel)
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
    #print "setting attributes for ",name
    if "a15" in name or "a20" in name or "a25" in name or "a30" in name or "a35" in name or "a40" in name:
        h = assignColor(h)
        hS.Add(h)
    if "ZL" in name or "ZTT" in name or "TTT" in name or "TTL" in name or "ZTL" in name or "jetFakes" in name:
        h = assignColor(h)
        hB.Add(h)




    return h,hS,hB,hD








if __name__ == "__main__":

    #Change these to input category files ?? via parser?
    from utils.Parametrization import Category
    from utils.Parametrization import Process
    import argparse

    parser = argparse.ArgumentParser(description="make full plots from root files containing histograms")
    #parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
    parser.add_argument("-i",  "--input", default="",  help="postfix string from previous MakeDataCard step")
    parser.add_argument("-o",  "--output", default="",  help="postfix string")
    parser.add_argument("-ch",  "--channel", default="mmmt",  help="postfix string")
    parser.add_argument("-c",  "--categories", default="categories.yaml",  help="categories yaml file")
    parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
    parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
    parser.add_argument("--dist", default="irBkg",  help="single distribution to plot")
    parser.add_argument("-mc",  "--mc", default=False,action='store_true',  help="Use only mc skip data")
    parser.add_argument("-mhs",  "--mhs", default=False,action='store_true',  help="make file containing histograms for datacards")
    parser.add_argument("-fh",  "--fh", default=False,action='store_true',  help="Make Finalized histograms")
    parser.add_argument("-ss",  "--signalScale", default=1.0,  help="Scale the Signal")
    args = parser.parse_args()

    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)




    #importing analysis categories and conventions
    with io.open(args.categories,'r') as catsyam:
        categories = yaml.load(catsyam)

    #loading fake factor and data driven methods
    with io.open(args.processes,'r') as prosyam:
        processes_special = yaml.load(prosyam)

    allcats={}

    for category in categories:
        #print category
        #print categories[category]['name']
        tempcat = Category()
        tempcat.name=categories[category]['name']
        tempcat.cuts=categories[category]['cuts']
        tempcat.newvariables=categories[category]['newvariables']
        tempcat.vars=categories[category]['vars']
        allcats[tempcat.name]=tempcat

    print "the categories   ",allcats
    newvars=[]
    variabledic={}

    #for newvar in  HAA_Inc_mmmt.newvariables.keys():
    for newvar in  allcats[args.channel+"_inclusive"].newvariables:
        newvars.append([newvar])
        variabledic[newvar]=[newvar,allcats[args.channel+"_inclusive"].newvariables[newvar][1],allcats[args.channel+"_inclusive"].newvariables[newvar][3],allcats[args.channel+"_inclusive"].newvariables[newvar][4]]
    variables = []
    for varHandle in allcats[args.channel+"_inclusive"].vars.keys():
        variables.append([varHandle])
        variabledic[varHandle]=allcats[args.channel+"_inclusive"].vars[varHandle]
    #variables = variables + newvars
    #print variables

    #The skimmed root file containing all the TTrees
    #file = ROOT.TFile(args.input,"read")
    fin = uproot.open(args.input)


    histolist = {}
    finalhists = {}
    processes = []
    histodict = {}
    dists = fin.keys()
    cats = [args.channel+"_inclusive"]
    vars = allcats[cats[0]].vars.keys()
    passvars = []

    for newvar in allcats[cats[0]].newvariables.keys():
        vars.append(newvar)

    #for ivar,var in enumerate(cat.vars.keys()):
    #for dist in fin.GetListOfKeys():
    systematics = ["Nominal","scale_m_etalt1p2Up"]
    for sys in systematics:
        histodict[sys]={}
        for distLong in dists:
            dist = distLong.split(";")[0]
            #distribution = dist.ReadObj()  #The TTree
            tree = fin[dist]
            masterArray = tree.arrays()
            histodict[sys][dist]={}
            for cat in cats:
                histodict[sys][dist][cat]={}
                for variableHandle in vars:
                    print variableHandle
                    if variableHandle in allcats[cats[0]].newvariables.keys():
                        var = variableHandle
                        bins = allcats[cat].newvariables[variableHandle][1]
                    else:
                        var = allcats[cat].vars[variableHandle][0]
                        bins = allcats[cat].vars[variableHandle][1]
                    if type(bins[0])==list:
                        histodict[sys][dist][cat][variableHandle] = ROOT.TH1D(str(dist),str(dist),bins[0][0],bins[0][1],bins[0][2])
                        try:
                            val = masterArray[var]
                            root_numpy.fill_hist(histodict[sys][dist][cat][variableHandle],val,masterArray["finalweight"])
                            passvars.append(var)
                        except:
                            print "problem with variable so skipping ",var
                            continue
                    else:
                        tmpbin = np.asarray(bins)
                        histodict[sys][dist][cat][variableHandle] = ROOT.TH1D(str(dist),str(dist),len(tmpbin)-1,tmpbin)
                        try:
                            val = masterArray[var]
                            root_numpy.fill_hist(histodict[sys][dist][cat][variableHandle],val,masterArray["finalweight"])
                            passvars.append(var)
                        except:
                            print "problem with variable so skipping ",var
                            continue
                    #histodict[sys][variablehandle+":"+allcats[cat].name+":"+dist] = root.TH1D(str(process),str(process),bins[0][0],bins[0][1],bins[0][2])
                    #root_numpy.fill_hist(histodict[variablehandle+":"+allcats[cat].name+":"+dist],val,masterArray["finalweight"])



    #for varCatDist in histodict.keys():

    #var is now the variable handle
    for sys in systematics:
        try:
            os.mkdir("outplots_"+args.output+"_"+sys)
        except:
            print "dir prob exists"

        for var in passvars:
            if args.mc:
                if var=="AMass":
                    fileout = open("outplots_"+args.output+"_"+sys+"/"+str(allcats[cats[0]].name)+"_info.txt","w")
                    fileout.write("Working on category "+allcats[cats[0]].name+"\n")
            else:
                if var=="AMass_blinded":
                    fileout = open("outplots_"+args.output+"_"+sys+"/"+str(allcats[cats[0]].name)+"_info.txt","w")
                    fileout.write("Working on category "+allcats[cats[0]].name+"\n")
                #if var in ["AMass","mll_m15","mll_m20","mll_m25","mll_m30","mll_m35","mll_m40","mll_m45","mll_m50","mll_m55","mll_m60"]:
                #    continue


            for cat in cats:
                #signal
                hSignals={}

                hDataDict={}
                hMCDict={}
                hMCFake1Dict={}
                hMCFake2Dict={}

                #print "divising MC into categories "
                hirBackground = ROOT.TH1F()

                hirBackground = histodict[sys][sys+"_"+args.dist][cat][var].Clone()


                #data

                histolist=[]

                hirBackground.SetLineColor(1)
                hirBackground.SetFillStyle(1001)
                hirBackground.SetFillColor(ROOT.TColor.GetColor("#13E2FE"))

                hirBackground.SetTitle("ZZ or ZH to 4l")

                hBkgTot = ROOT.THStack()
                hBkgTot.Add(hirBackground)










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

                doRatio=False

                if not doRatio:
                    c.SetLeftMargin(L/W)
                    c.SetRightMargin(R/W)
                    c.SetTopMargin(T/H)
                    c.SetBottomMargin(B/H)
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
                lumi=add_lumi("2016")
                cms=add_CMS()
                pre=add_Preliminary(args.channel)
                xR=0.65
                l=ROOT.TLegend(xR,0.55,xR+0.28,0.9);
                if var=="mll_fine":
                    l=ROOT.TLegend(0.40,0.55,0.40+0.28,0.9);
                l.AddEntry(hirBackground)

                #print "signal entries   ",hSignal.GetEntries()
                hirBackground.GetXaxis().SetTitle(variabledic[var][3]+variabledic[var][2])
                hirBackground.GetYaxis().SetTitle("Events")
                hirBackground.SetTitle("")
                hirBackground.Draw("HIST")

                #hirBackground.Draw("same")
                #for hist in histolist:
                #    hist.Draw("same")
                lumi.Draw()
                cms.Draw()
                pre.Draw(args.channel)
                l.Draw()

                #TPad 2 for ratio
                c.cd()

                #print "with cuts ",allcats[cati].cuts
                #print "data entries ",hData.GetEntries()
                #print "background entries ",hBackground.GetEntries()

                c.SaveAs("outplots_"+args.output+"_"+sys+"/"+var+"_"+str(cats[0])+".png")

                hirBackground.Delete()
