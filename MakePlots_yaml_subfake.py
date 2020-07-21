#!/usr/bin/env python
import ROOT
import re
import math
from array import array
#from collections import OrderedDict
#import varCfgPlotter
import argparse
import os
import io
import yaml
#from HttStyles import GetStyleHtt
#from HttStyles import MakeCanvas

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
    from utils.Parametrization import Category
    from utils.Parametrization import Process
    import argparse

    parser = argparse.ArgumentParser(description="make full plots from root files containing histograms")
    #parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
    parser.add_argument("-i",  "--input", default="",  help="postfix string from previous MakeDataCard step")
    parser.add_argument("-o",  "--output", default="",  help="postfix string")
    parser.add_argument("-c",  "--categories", default="categories.yaml",  help="categories yaml file")
    parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
    parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
    parser.add_argument("-dd",  "--datadriven", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-mc",  "--mc", default=False,action='store_true',  help="Use only mc skip data")
    args = parser.parse_args()

    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)



    #Gather the analysis datasets and info 
    sampleDict = {}
 
    #with open("MCsamples_2016_v6_yaml.csv")  as csvfile:
    #    reader = csv.reader(csvfile, delimiter=',')
    #    for row in reader:
    #        sampleDict[row[0]] = [row[1],row[2],row[3],row[4],row[5],row[6]]

    for line in open(args.csvfile,'r').readlines() :
            #[nickname]        = [category,xsec,numberOfEvents,finishedEvents,idk?,DASDataset]
            if len(line.split(',')[2].split("*"))>1:
                tempval=1.0
                for val in line.split(',')[2].split("*"):
                    tempval = tempval * float(val)
                row = line.split(',')
                sampleDict[row[0]] = [row[1],tempval,row[3],row[4],row[5],row[6]]
            else:
                row = line.split(',')
                sampleDict[row[0]] = [row[1],float(row[2]),row[3],row[4],row[5],row[6]]

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

    newvars=[]
    variabledic={}

    #for newvar in  HAA_Inc_mmmt.newvariables.keys():
    for newvar in  allcats[allcats.keys()[0]].newvariables:
        newvars.append([newvar])
        variabledic[newvar]=[newvar,allcats[allcats.keys()[0]].newvariables[newvar][1],allcats[allcats.keys()[0]].newvariables[newvar][3],allcats[allcats.keys()[0]].newvariables[newvar][4]]
    variables = [] 
    for varHandle in allcats[allcats.keys()[0]].vars.keys():
        variables.append([varHandle])
        variabledic[varHandle]=allcats[allcats.keys()[0]].vars[varHandle]
    variables = variables + newvars
    #print variables
    filelist = {}
    processes = []
    
    for var in variables:
        filelist[var[0]] = var[0]+".root"



    #loading standard processes
    HAA_processes={}
    for sample in sampleDict.keys():
        #print processes[process]
        temppro = Process() 
        temppro.nickname=sample
        temppro.file=sample+"_2016.root"
        temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
        temppro.cuts={sampleDict[sample][0]:""}
        if "ggTo2mu2tau" in sample:
            temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(137.5*31.05*0.00005)}
        if "W" in sample and "Jets" in sample: 
            temppro.file=sample+"_2016.root"
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
            temppro.cuts={"W":"","WL":[["gen_match_4",">=",5]],"WJ":[["gen_match_4",">",5]]}
        if "TT" in sample and not "TTTT" in sample and not "TTHH" in sample: 
            temppro.file=sample+"_2016.root"
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3]}
            temppro.cuts={"TT":"","TTT":[["gen_match_4","==",5]],"TTL":[["gen_match_4",">=",5]],"TTJ":[["gen_match_4",">",5]]}
        HAA_processes[temppro.nickname]=temppro

    #loading special processes ... fake factor and data
    for process in processes_special:
        temppro = Process() 
        temppro.nickname=processes_special[process]['nickname']
        temppro.cuts=processes_special[process]['cuts']
        temppro.weights=processes_special[process]['weights']
        temppro.file=processes_special[process]['file']
        HAA_processes[temppro.nickname]=temppro

    for pro in HAA_processes.values():
        for subpro in pro.cuts.keys():
            processes.append(subpro)
    processes=list(dict.fromkeys(processes))
    print "Will plot the following physics processes ",processes


    for ivar,var in enumerate(filelist.keys()):
    #for ivar,var in enumerate(cat.vars.keys()):
        print "Working on file ",filelist[var],"  for variable   ",var
        try:
            file = ROOT.TFile("out"+str(args.input)+"/"+filelist[var],"read")
        except:
            print "no file for ",filelist[var]
            continue
        #bins=HAA_Inc_mmmt.binning + HAA
         
        #for catname in enumerate():
        cati = 0
        #for categoryKey in file.GetListOfKeys():
        for categoryKey in file.GetListOfKeys():
            category = categoryKey.ReadObj()
            print "Working on category ",category.GetName()
            try:
                os.mkdir("outplots_"+args.output+"_"+str(category.GetName()))
            except:
                print "dir prob exists"
            
            if args.mc:
                if var=="AMass":
                    fileout = open("outplots_"+args.output+"_"+str(category.GetName())+"/"+str(allcats[category.GetName()].name)+"_info.txt","w")  
                    fileout.write("Working on category "+allcats[category.GetName()].name+"\n")
            else: 
                if var=="AMass_blinded":
                    fileout = open("outplots_"+args.output+"_"+str(category.GetName())+"/"+str(allcats[category.GetName()].name)+"_info.txt","w")  
                    fileout.write("Working on category "+allcats[category.GetName()].name+"\n")
                #if var in ["AMass","mll_m15","mll_m20","mll_m25","mll_m30","mll_m35","mll_m40","mll_m45","mll_m50","mll_m55","mll_m60"]:
                #    continue

            #signal 
            hSignal = category.Get("a40")
            
            hDataDict={}
            hMCDict={}
            hMCFake1Dict={}
            hMCFake2Dict={}

            #obtaining the histograms
            for histokey in category.GetListOfKeys():
                histo = histokey.ReadObj()
                histoname = histo.GetName()
                print histoname.split("fake1_")[0]
                if histoname not in ["data_obs","FF_1","FF_2","FF_12"] and histoname.split("fake1_")[0]=="":
                    hMCFake1Dict[histoname]=histo
                if histoname not in ["data_obs","FF_1","FF_2","FF_12"] and histoname.split("fake2_")[0]=="":
                    hMCFake2Dict[histoname]=histo
                #if histoname not in ["data_obs","FF_1","FF_2","FF_12","TTL","WL","ZL","WJ","TTJ","TTT","FF","prompt1","prompt2"] and histoname.split("fake1_")[0]!="" and histoname.split("fake2_")[0]!="":
                if histoname in ["DY","W","TT","ST","EWK","ZZ","ZHToTauTau","vbf","WHTT","ttZ","ttW","WWZ","WZZ","ZZZ","WWW_4F","HZJ","Other","rare","WZ"] and histoname.split("fake1_")[0]!="" and histoname.split("fake2_")[0]!="":
                    hMCDict[histoname]=histo
                elif histoname in ["data_obs","FF_1","FF_2","FF_12"]:
                    hDataDict[histoname]=histo
            print hMCDict
            #print hMCFake1Dict

            #subtracting the fakes 
            print "subtracting fakes "
            for histo in hMCDict.keys():
                hMCDict[histo].Add(hMCFake1Dict["fake1_"+str(histo)],-1)
                hMCDict[histo].Add(hMCFake2Dict["fake2_"+str(histo)],-1)

            print "divising MC into categories "
            hBackground = ROOT.TH1F()
            hirBackground = ROOT.TH1F()
            hrareBackground = ROOT.TH1F()
            h3alphBackground = ROOT.TH1F()

            Bkg = ["DY","W","TT","ST","EWK"]
            irBkg = ["ZZ","ZHToTauTau","vbf","WHTT"]
            TrialphaBkg = ["ttZ","ttW","WWZ","WZZ","ZZZ","WWW_4F","HZJ"]
            rareBkg = ["Other","rare","WZ"]

            for histo in hMCDict.keys():
                if histo in Bkg:
                    if hBackground.GetName()=="":
                        hBackground = hMCDict[histo].Clone()
                    else:
                        hBackground.Add(hMCDict[histo])
                if histo in irBkg:
                    if hirBackground.GetName()=="":
                        hirBackground = hMCDict[histo].Clone()
                    else:
                        hirBackground.Add(hMCDict[histo])
                if histo in TrialphaBkg:
                    if h3alphBackground.GetName()=="":
                        h3alphBackground = hMCDict[histo].Clone()
                    else:
                        h3alphBackground.Add(hMCDict[histo])
                if histo in rareBkg:
                    if hrareBackground.GetName()=="":
                        hrareBackground = hMCDict[histo].Clone()
                    else:
                        hrareBackground.Add(hMCDict[histo])

            #hrareBackground.Add(category.Get("GluGluTo2mu2nu"))
            #hrareBackground.Add(category.Get("GluGluTo4mu"))

            #data
            if not args.mc:
                hData = category.Get("data_obs")
                dataMax = hData.GetMaximum()*1.25
                dataMin = hData.GetMinimum()*1.25
                if var in ["mll_fine"]:
                    hData.SetBinContent(7,-999)
                    hData.SetBinContent(8,-999)
                    hData.SetBinContent(9,-999)
                    hData.SetBinContent(10,-999)
                    hData.SetBinContent(11,-999)
                if var in ["mll"]:
                    hData.SetBinContent(3,-999)
                    hData.SetBinContent(4,-999)
                    hData.SetBinContent(5,-999)
                    hData.SetBinContent(6,-999)

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
            hBackground.SetTitle("DY+Jets, T#bar{T}, ST , EWK")
            #hirBackground.SetTitle("irreducible")
            hirBackground.SetTitle("ZZ or ZH to 4l")
            h3alphBackground.SetTitle("3 #alpha")
            hrareBackground.SetTitle("TTXX, WZ to Inv, 4#alpha")

            hBkgTot = ROOT.THStack()
            hBkgTot.Add(hirBackground)
            hBkgTot.Add(hrareBackground)
            hBkgTot.Add(h3alphBackground)

            if args.datadriven:
                hFF1 = category.Get("FF_1")
                hFF2 = category.Get("FF_2")
                hFF12 = category.Get("FF_12")
                hFF = hFF1.Clone()
                #hFF.SetTitle("jetFakes")
                hFF.SetTitle("Jet faking #tau")
                hFF.Add(hFF2)
                hFF.Add(hFF12,-1)
                hFF.Add(hFF12,-1)
                #hFF.Scale(0.5)
                hFF.SetLineColor(1)
                hFF.SetFillStyle(1001)
                hFF.SetFillColor(ROOT.TColor.GetColor("#CF8AC8"))

            if args.datadrivenZH:
                hData = category.Get("data_obs")
                hPrompt = category.Get("prompt")
                #hFake1 = category.Get("fake1")
                #hFake2 = category.Get("fake2")
                #hirBackground.Add(hFake1,-1)
                #hirBackground.Add(hFake2,-1)
                hFF1 = category.Get("FF_1")
                hFF2 = category.Get("FF_2")
                hFF12 = category.Get("FF_12")
                hFF = hFF1.Clone()
                #hFF.SetTitle("jetFakes")
                hFF.SetTitle("Jet faking #tau")
                hFF.Add(hFF2)
                hFF.Add(hFF12,-1)
                #hFF.Add(hPrompt,-1)
                #hFF.Add(hFF12,-1)
                #hFF.Scale(0.5)
                hFF.SetLineColor(1)
                hFF.SetFillStyle(1001)
                hFF.SetFillColor(ROOT.TColor.GetColor("#CF8AC8"))

            

            if not (args.datadriven or args.datadrivenZH):
                hbkg = hBackground.Clone()
                hBkgTot.Add(hBackground)
            else:
                hbkg = hFF.Clone()
                hBkgTot.Add(hFF)

            
            hbkg.Add(hirBackground)
            hbkg.Add(hrareBackground)
            hbkg.Add(h3alphBackground)
                



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
            pre=add_Preliminary()
            xR=0.65
            l=ROOT.TLegend(xR,0.55,xR+0.28,0.9);
            if var=="mll_fine":
                l=ROOT.TLegend(0.40,0.55,0.40+0.28,0.9);
            if not args.mc:
                l.AddEntry(hData)
            l.AddEntry(hSignal)
            if not (args.datadriven or args.datadrivenZH):
                l.AddEntry(hBackground)
            else:
                l.AddEntry(hFF)
            l.AddEntry(h3alphBackground)
            l.AddEntry(hirBackground)
            l.AddEntry(hrareBackground)

            #print "signal entries   ",hSignal.GetEntries()
            hSignal.GetXaxis().SetTitle(variabledic[var][3]+variabledic[var][2])
            hSignal.GetYaxis().SetTitle("Events")
            #hBkgTot.GetXaxis().SetTitle(variabledic[var][3]+variabledic[var][2])
            #hBkgTot.GetYaxis().SetTitle("Events")
            if not args.mc:
                hData.SetTitle("")
                hData.SetMarkerStyle(8)
                hData.GetXaxis().SetTitle(variabledic[var][3]+variabledic[var][2])
                hData.GetYaxis().SetTitle("Events")
                hData.GetYaxis().SetRangeUser(dataMin,dataMax)
                hData.Draw("ep")
                hBkgTot.Draw("HIST same")
                hSignal.Draw("HIST same")
                hData.Draw("ep same")
            else:
                if(hBkgTot.GetMaximum() > hSignal.GetMaximum()):
                    hBkgTot.Draw("HIST")
                    hSignal.Draw("HIST same")
                else:
                    hSignal.SetTitle("")
                    hSignal.Draw("HIST")
                    hBkgTot.Draw("HIST same")
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
            if not args.mc:
                hData2=hData.Clone()
                pad2.Draw()
                pad2.cd()
                hData2.Divide(hbkg)
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
            else:
                hSignal2=hSignal.Clone()
                pad2.Draw()
                pad2.cd()
                hSignal2.Divide(hbkg)
                hSignal2.SetMarkerStyle(20)
                hSignal2.SetTitleSize  (0.12,"Y")
                hSignal2.SetTitleOffset(0.40,"Y")
                hSignal2.SetTitleSize  (0.12,"X")
                hSignal2.SetLabelSize  (0.10,"X")
                hSignal2.SetLabelSize  (0.08,"Y")
                #hSignal2.GetYaxis().SetRangeUser(0.62,1.38)
                hSignal2.GetYaxis().SetRangeUser(0.3,2.2)
                hSignal2.GetYaxis().SetNdivisions(305)
                hSignal2.GetYaxis().SetTitle("Sig/Bkg   ")
                hSignal2.Draw("P");
                pad2.Draw()

            #print "with cuts ",allcats[cati].cuts
            #print "data entries ",hData.GetEntries()
            #print "background entries ",hBackground.GetEntries()
            if args.mc:
                if var=="AMass":
                    for key in allcats[category.GetName()].cuts.keys():
                        for cut in allcats[category.GetName()].cuts[key]:
                            fileout.write(str(cut))
                        fileout.write("\n")
                    fileout.write("signal entries "+str(hSignal.Integral())+"\n")
                    fileout.write("background entries "+str(hBackground.Integral())+"\n")
            else:
                if var=="AMass_blinded":
                    for key in allcats[category.GetName()].cuts.keys():
                        for cut in allcats[category.GetName()].cuts[key]:
                            fileout.write(str(cut))
                        fileout.write("\n")
                    fileout.write("data entries "+str(hData.Integral())+"\n")
                    fileout.write("background entries "+str(hbkg.Integral())+"\n")
            
            c.SaveAs("outplots_"+args.output+"_"+str(category.GetName())+"/"+var+"_"+str(category.GetName())+".png")
            cati = cati +1

            hBackground.Delete()
            hirBackground.Delete()
            h3alphBackground.Delete()
            if not args.mc:
                hData.Delete()
                hData2.Delete()
            else:
                hSignal2.Delete()
            hSignal.Delete()
            hBkgTot.Delete()
            hbkg.Delete()


