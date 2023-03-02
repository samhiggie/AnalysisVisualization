# !/usr/bin/env python
#########################
#Author: Sam Higginbotham
'''

* File Name : CreateFakeHistos_v7.py

* Purpose : creating ratio files to measure probability of objects faking other objects... important for tau physics ;)

* Usage:
python CreateFakeHistos_v7.py -c cat_ff_Z_e.yaml -csv MCsamples_2018_paper.csv -i /eos/home-s/shigginb/HAA_ntuples/legacy_2018/ -p processes_loose_mmet.yaml -o 2018_ff_Z_e -ch mmet -ff -fi 2018_ff_Z_dm_mmet -fo 2018_ff_Z_mmet -year 2018

* Creation Date : Dec-07-2022


'''
#########################
import sys
import os
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kFatal
import uproot
import pandas
import root_numpy
import numpy as np
from array import array
import itertools
import operator
import csv
import datetime
import time
import threading
import yaml
import glob
#user definitions
from utils.Parametrization import *


#global vars for setting up the analysis categories
allcats={}
HAA_processes={}
finalDistributions={}
filelist = {}
weightHistoDict={}
jetWeightMultiplicity={}
EventWeights={}
datadrivenPackage={}

from MakeDistributions import initialize
import argparse


parser = argparse.ArgumentParser(description="This file generates root files containing Histograms ... files in utils contain selections and settings")
parser.add_argument("-o",  "--outname", default="",  help="postfix string")
parser.add_argument("-fi",  "--ffin", default="",  help="fake factor files")
parser.add_argument("-year",  "--year", default="2016",  help="Year")
parser.add_argument("-fo",  "--ffout", default="",  help="fake factor files to output")
parser.add_argument("-ft",  "--ftype", default="",  help="")
parser.add_argument("-lt",  "--lepton", default="",  help="")
parser.add_argument("-c",  "--categories", default="categories_array.yaml",  help="categories yaml file")
parser.add_argument("-catt",  "--cattight", default="",  help="if measuring single leg pick one tight id/iso category")
parser.add_argument("-catl",  "--catloose", default="",  help="if measuring single leg pick one loose id/iso category")
parser.add_argument("-ch",  "--channel", default="mmmt",  help="Please list the channel for fake factor histograms")
parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
parser.add_argument("-ddHAA",  "--datadrivenHAA", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-i",  "--dir", default="/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov7_basic_10_6_4/src/2016_v7/",  help="Input files")
parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
parser.add_argument("-pt",  "--prompt", default=False,action='store_true',  help="When making fake rate histograms, estimate the prompt contribution")
parser.add_argument("-spt",  "--subprompt", default=False,action='store_true',  help="When making fake rate histograms, subtract the prompt contribution")
parser.add_argument("-dm",  "--datameasure", default=False,action='store_true',  help="Use DataDriven Method measure part")
parser.add_argument("-dmZH",  "--datameasureZH", default=False,action='store_true',  help="Use DataDriven Method measure part")
parser.add_argument("-dmOD",  "--onlydata", default=False,action='store_true',  help="only produced distributions derived from data ntuple.")
parser.add_argument("-ex",  "--extract", default=False,action='store_true',  help="Additional Cuts for Extraction")
parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-ddSM",  "--datadrivenSM", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-sys",  "--systematics", default=False,action='store_true',  help="Run on systematic trees, don't run nominal at the same time")
parser.add_argument("-ff",  "--makeFakeHistos", default=False,action='store_true',  help="Just make fake rate histos")

args = parser.parse_args()
print(args)

def add_lumi(lowX,lowY,year,doRatio=True):
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.04)
    lumi.SetTextFont (   42 )
    if year == "2016":
        lumi.AddText(str(year)+" 35.9 fb^{-1} (13 TeV)")
    if year == "2017":
        lumi.AddText(str(year)+" 41.8 fb^{-1} (13 TeV)")
    if year == "2018":
        lumi.AddText(str(year)+" 59.7 fb^{-1} (13 TeV)")
    if year == "RunII":
        lumi.AddText(str(year)+" 137 fb^{-1} (13 TeV)")
    return lumi

def add_CMS(lowX,lowY,doRatio=True):
    print("+++++++++++++++++++++++++++++++++++++++++++")
    print(lowX,lowY)
    print("+++++++++++++++++++++++++++++++++++++++++++")
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.06)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS")
    return lumi

def add_Preliminary(lowX,lowY,channel="mmmt", doRatio=True):

    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(52)
    lumi.SetTextSize(0.04)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    #lumi.AddText("Preliminary "+channel)
    lumi.AddText("Preliminary ")
    return lumi
def setUpPlot(ratio):
    H = 600
    W = 600
    
    H_ref = 600
    W_ref = 600


    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.16*W_ref
    R = 0.04*W_ref


    xR=0.65
    B_ratio = 0.1*H_ref
    T_ratio = 0.03*H_ref

    B_ratio_label = 0.3*H_ref

    c=ROOT.TCanvas("canvas","",0,0,600,600)
    pad1 = ROOT.TPad("pad1","pad1",0.0016,0.291,1.0,1.0)
    pads = []

    if not ratio:
        c.SetLeftMargin(L/W)
        c.SetRightMargin(R/W)
        c.SetTopMargin(T/H)
        c.SetBottomMargin(B/H)
        pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
        L = 0.13*W_ref
        pad1.SetLeftMargin     (L/W)
        pad1.SetRightMargin    (R/W)
        pad1.SetTopMargin      (T/H)
        pad1.SetBottomMargin   (B/H)

        l=ROOT.TLegend(xR,0.45,xR+0.25,0.9);
        pads.append(pad1)
    else:
        L = 0.16*W_ref
        pad1.SetTicks(0,0)
        pad1.SetLeftMargin(L/W)
        pad1.SetRightMargin(R/W)
        pad1.SetTopMargin(T/H)
        pad1.SetBottomMargin(B_ratio/H)
        pad1.SetFillColor(0)
        pad1.SetBottomMargin(0)

        pad2 = ROOT.TPad("pad2","The lower pad",0,0,1.0,0.29)
        pad2.SetLeftMargin(L/W)
        pad2.SetRightMargin(R/W)
        pad2.SetTopMargin(T_ratio/H)
        pad2.SetTopMargin(0.007)
        pad2.SetBottomMargin(B_ratio_label/H)
        pad2.SetGridy(1)
        pad2.SetFillColor(4000)
        l=ROOT.TLegend(xR,0.55,xR+0.28,0.9);
        pads.append(pad1,pad2)

    return c,pads,l

def getPlotStuffs(year,channel):
    doRatio=False
    Ylow=0.735
    Xlow=0.58
    lumi=add_lumi(Xlow,Ylow,year,doRatio)
    Xlow=0.27
    cms=add_CMS(Xlow,Ylow,doRatio)
    Xlow=0.39
    pre=add_Preliminary(Xlow,Ylow,channel, doRatio)
    return lumi,cms,pre
def createFakeFactorHistoslocal(allcats, inputFFile,prompt,ftype):
    newVarVals={}
    if prompt:
        treetypes = ["Nominal_data_obs","Nominal_prompt1","Nominal_prompt2"]
    else:
        treetypes = ["Nominal_data_obs"]
    histodict = {c:{t: dict() for t in treetypes} for c in allcats.keys()}
    print("categories  ",histodict.keys())


    #creating the Fake Factor histograms from pre-defined numpy arrays
    for cat in histodict.keys():
        #print cat
        if "inclusive" in cat: continue
        for treename in treetypes:
            tree = inputFFile[cat][treename]
            fakefactorArray = tree.arrays(library="np")
            #print fakefactorArray.keys()
            if "prompt" in treename: 
                if ftype=="firstleg":
                    mask = fakefactorArray["gen_match_3"] == 5
                if ftype=="secondleg":
                    mask = fakefactorArray["gen_match_4"] == 5
                #print fakefactorArray["gen_match_3"]
                #print mask
                for var,branch in fakefactorArray.items():
                    fakefactorArray[var] = branch[mask]
                
            #print("gen_match",fakefactorArray["gen_match_3"])
            #print "cats ",cat,"  tree name ",treename
            #print "entries ",len(fakefactorArray["finalweight"])
            #vars = fakefactorArray.keys()
            vars = allcats[cat].vars
            for variableHandle in vars:
                #print variableHandle
                var = allcats[cat].vars[variableHandle][0]
                #var = variableHandle
                #print var,"  ",variableHandle
                if "jpt" in variableHandle and "ff" in variableHandle:
                    firstnum = variableHandle.split("_")[1]
                    secondnum = str(int(firstnum)-1)
                    var = variableHandle+"_"+secondnum
                    #print var
                bins = allcats[cat].vars[variableHandle][1]
                if type(bins[0])==list:
                    histodict[cat][treename][variableHandle] = ROOT.TH1D(str(variableHandle),str(variableHandle),bins[0][0],bins[0][1],bins[0][2])
                    val = fakefactorArray[var]
                    #val = returnArray(fakefactorArray,var)
                    root_numpy.fill_hist(histodict[cat][treename][variableHandle],val,fakefactorArray["finalweight"])
                else:
                    tmpbin = np.asarray(bins)
                    histodict[cat][treename][variableHandle] = ROOT.TH1D(str(variableHandle),str(variableHandle),len(tmpbin)-1,tmpbin)
                    val = fakefactorArray[var]
                    #val = returnArray(fakefactorArray,var)
                    root_numpy.fill_hist(histodict[cat][treename][variableHandle],val,fakefactorArray["finalweight"])
    #now I have the histograms ... so using standard methods
    return histodict

if __name__ == "__main__":
    import datetime
    begin_time = datetime.datetime.now()
    from multiprocessing import *
    import multiprocessing as mp
    import logging
    import os
    import shutil


    #trial of global settings
    #import AnalysisSetup as AS
    from utils.Weights import jet_exc_samples



    print("setting up all the expect distributions")

    allcats,HAA_processes,finalDistributions,\
    weightHistoDict,jetWeightMultiplicity,datadrivenPackage,fakefactorObj = initialize(args)

    

    print("MUST HAVE CREATED FF DISTRIBUTIONS BEFORE!")



    try:
        os.mkdir("FFhistos_"+str(args.ffout))
    except:
        pass

    histodict = {}
    #inputFFile = ROOT.TFile.Open("FFskim_"+str(args.ffin)+".root","read")
    infofile = open("FFhistos_"+str(args.ffout)+"/fake_"+str(args.ftype)+".txt","w")

    with uproot.open("skimmed_"+str(args.ffin)+".root") as inputFFile:
        for cat in inputFFile.keys():
            infofile.write(cat)
            infofile.write("\n")
        histodict = createFakeFactorHistoslocal(allcats, inputFFile,args.prompt,args.ftype)

    datadrivenPackage={}
    datadrivenPackage["bool"]=args.datadrivenZH
    fakemeasurefile = ROOT.TFile.Open("FFhistos_"+str(args.ffout)+"/fake_"+str(args.ftype)+".root","RECREATE")
    fakemeasurefile.cd()

    #another try end of year 2022
    #ss_1_tight = histodict[args.channel+"_FF_SS_validation"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_loose = histodict[args.channel+"_FF_SS_1_loose"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_2_tight = histodict[args.channel+"_FF_SS_validation"]["Nominal_data_obs"]["pt_4_ff"]
    #ss_2_loose = histodict[args.channel+"_FF_SS_2_loose"]["Nominal_data_obs"]["pt_4_ff"]
    #validation too small yield... loosen other leg ??
    #ss_1_tight = histodict[args.channel+"_FF_SS_2_loose"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_loose = histodict[args.channel+"_FF_SS_12_loose"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_2_tight = histodict[args.channel+"_FF_SS_1_loose"]["Nominal_data_obs"]["pt_4_ff"]
    #ss_2_loose = histodict[args.channel+"_FF_SS_12_loose"]["Nominal_data_obs"]["pt_4_ff"]

    #inclusive definitions of cuts
    #ss_1_tight = histodict[args.channel+"_FF_SS_2_loose_inc"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_loose = histodict[args.channel+"_FF_SS_12_loose_inc"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_2_tight = histodict[args.channel+"_FF_SS_1_loose_inc"]["Nominal_data_obs"]["pt_4_ff"]
    #ss_2_loose = histodict[args.channel+"_FF_SS_12_loose_inc"]["Nominal_data_obs"]["pt_4_ff"]
    #ss_1_tight = histodict[args.channel+"_FF_SS_validation"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_loose = histodict[args.channel+"_FF_SS_1_loose_inc"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_2_tight = histodict[args.channel+"_FF_SS_validation"]["Nominal_data_obs"]["pt_4_ff"]
    #ss_2_loose = histodict[args.channel+"_FF_SS_2_loose_inc"]["Nominal_data_obs"]["pt_4_ff"]



    #haa old
    #just electon mmet 
    #ss_1_tight = histodict[args.channel+"_FF_SS_2_loose"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_loose = histodict[args.channel+"_FF_SS_12_loose_inc"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_tight_prompt = histodict[args.channel+"_FF_SS_2_loose"]["Nominal_prompt1"]["pt_3_ff"]
    #ss_1_loose_prompt = histodict[args.channel+"_FF_SS_12_loose_inc"]["Nominal_prompt1"]["pt_3_ff"]
    #ss_1_tight = histodict[args.channel+"_FF_SS_2_loose"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_loose = histodict[args.channel+"_FF_SS_12_loose"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_tight_prompt = histodict[args.channel+"_FF_SS_2_loose"]["Nominal_prompt1"]["pt_3_ff"]
    #ss_1_loose_prompt = histodict[args.channel+"_FF_SS_12_loose"]["Nominal_prompt1"]["pt_3_ff"]


    #ss_1_tight_prompt = histodict[args.channel+"_FF_SS_validation"]["Nominal_prompt1"]["pt_3_ff"]
    #ss_1_loose_prompt = histodict[args.channel+"_FF_SS_1_loose_inc"]["Nominal_prompt1"]["pt_3_ff"]

    #ss_2_tight = histodict[args.channel+"_FF_SS_1_loose"]["Nominal_data_obs"]["pt_4_ff"]
    #ss_2_loose = histodict[args.channel+"_FF_SS_12_loose_inc"]["Nominal_data_obs"]["pt_4_ff"]

    #just muon mmmt
    #ss_1_tight = histodict[args.channel+"_FF_SS_2_loose_inc"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_loose = histodict["mmmt_FF_SS_1_haa_loose_inc"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_tight_prompt = histodict[args.channel+"_FF_SS_2_loose_inc"]["Nominal_prompt1"]["pt_3_ff"]
    #ss_1_loose_prompt = histodict["mmmt_FF_SS_1_haa_loose_inc"]["Nominal_prompt1"]["pt_3_ff"]
    #ss_1_tight = histodict[args.channel+"_FF_SS_2_loose_inc"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_loose = histodict["mmmt_FF_SS_1_haa_loose"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_tight_prompt = histodict[args.channel+"_FF_SS_2_loose_inc"]["Nominal_prompt1"]["pt_3_ff"]
    #ss_1_loose_prompt = histodict["mmmt_FF_SS_1_haa_loose"]["Nominal_prompt1"]["pt_3_ff"]
    #2 in use
    #ss_1_tight = histodict[args.channel+"_FF_SS_validation"]["Nominal_data_obs"]["pt_3_ff"]
    #ss_1_loose = histodict["mmmt_FF_SS_1_haa_loose"]["Nominal_data_obs"]["pt_3_ff"]

    #ss_1_tight_prompt = histodict[args.channel+"_FF_SS_validation"]["Nominal_prompt1"]["pt_3_ff"]
    #ss_1_loose_prompt = histodict["mmmt_FF_SS_1_haa_loose"]["Nominal_prompt1"]["pt_3_ff"]

    #ss_2_tight = histodict[args.channel+"_FF_SS_1_loose"]["Nominal_data_obs"]["pt_4_ff"]
    #ss_2_loose = histodict[args.channel+"_FF_SS_12_loose_inc"]["Nominal_data_obs"]["pt_4_ff"]



    #ss_1_tight.Add(ss_1_tight_prompt,-1)
    #ss_2_tight.Add(ss_2_tight_prompt,-1)
    #ss_1_loose.Add(ss_1_loose_prompt,-1)
    #ss_2_loose.Add(ss_2_loose_prompt,-1)
    if args.ftype=="firstleg" or args.ftype=="secondleg":
        if args.ftype=="firstleg":  legnum = str(3)
        if args.ftype=="secondleg": legnum = str(4)

        #cattight = "mmmt_FF_SS_validation"
        #cattight = "mmmt_FF_SS_1_tightnoiso"

        #cattight = "mmet_FF_SS_validation"
        cattight = args.cattight

        #cattight = "mmet_FF_SS_1_tightnoiso"
        #cattight = "mmet_FF_SS_1_tightiso"

        #catloose = "mmmt_FF_SS_1_haa_loose"
        #catloose = "mmmt_FF_SS_1_haa_loose"

        #catloose = "mmet_FF_SS_1_loose_inc"
        catloose = args.catloose

        #catloose = "mmet_FF_SS_1_looseiso"

        measurecats = [cattight,catloose]
        for cat in measurecats:
            infofile.write(cat)
            infofile.write("\n")
            #print(allcats[cat].cuts)
            for i,cut in allcats[cat].cuts.items():
                infofile.write(i)
                infofile.write(" : ")
                for subcut in cut:
                    infofile.write(str(subcut))
                infofile.write("\n")

        ss_1_tight = histodict[cattight]["Nominal_data_obs"]["pt_"+legnum+"_ff"]
        ss_1_loose = histodict[catloose]["Nominal_data_obs"]["pt_"+legnum+"_ff"]

        if args.prompt:
            ss_1_tight_prompt = histodict[cattight]["Nominal_prompt1"]["pt_"+legnum+"_ff"]
            ss_1_loose_prompt = histodict[catloose]["Nominal_prompt1"]["pt_"+legnum+"_ff"]

        print("data in tight ",ss_1_tight.Integral())
        infofile.write("data in tight \n")
        infofile.write(str(ss_1_tight.Integral()))
        infofile.write("\n")

        if args.prompt:
            print("prompt in tight ",ss_1_tight_prompt.Integral())
            infofile.write("prompt in tight "+str(ss_1_tight_prompt.Integral())+"\n")

        print("data in loose ",ss_1_loose.Integral())
        infofile.write("data in loose \n")
        infofile.write(str(ss_1_loose.Integral()))
        infofile.write("\n")

        if args.prompt:
            print("prompt in loose ",ss_1_loose_prompt.Integral())
            infofile.write("prompt in loose "+str(ss_1_loose_prompt.Integral())+"\n")

        if args.subprompt:
            ss_1_tight.Add(ss_1_tight_prompt,-1)
            ss_1_loose.Add(ss_1_loose_prompt,-1)
            print("data after subtraction in tight ",ss_1_tight.Integral())
            infofile.write("data after subtraction in tight \n")
            infofile.write(str(ss_1_tight.Integral()))
            infofile.write("\n")
            print("data after subtraction in loose ",ss_1_loose.Integral())
            infofile.write("data after subtraction in loose \n")
            infofile.write(str(ss_1_loose.Integral()))
            infofile.write("\n")

        c = ROOT.TCanvas()
        pad = ROOT.TPad()
        padlist=[]
        leg = ROOT.TLegend()

        c,padlist,leg = setUpPlot(False)
        pad = padlist[0]
        histodict[cattight]["Nominal_data_obs"]["pt_"+legnum+"_ff_fine"].SetTitle("Leg 1 Lepton p_T in Numerator")
        histodict[cattight]["Nominal_data_obs"]["pt_"+legnum+"_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[cattight]["Nominal_data_obs"]["pt_"+legnum+"_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[cattight]["Nominal_data_obs"]["pt_"+legnum+"_ff_fine"].Draw("PE")
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leg1_leptonPt_Numerator_"+args.channel+"_"+args.ftype+".png")
        c.SetLogy(1)
        histodict[cattight]["Nominal_data_obs"]["pt_"+legnum+"_ff_fine"].Draw("PE")
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leg1_leptonPt_Numerator_"+args.channel+"_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)

        histodict[catloose]["Nominal_data_obs"]["pt_"+legnum+"_ff_fine"].SetTitle("Leg 1 Lepton p_T in Denominator")
        histodict[catloose]["Nominal_data_obs"]["pt_"+legnum+"_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[catloose]["Nominal_data_obs"]["pt_"+legnum+"_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[catloose]["Nominal_data_obs"]["pt_"+legnum+"_ff_fine"].Draw("PE")
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leg1_leptonPt_Denominator_"+args.channel+"_"+args.ftype+".png")

        c.SetLogy(1)
        histodict[catloose]["Nominal_data_obs"]["pt_"+legnum+"_ff_fine"].Draw("PE")
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leg1_leptonPt_Denominator_"+args.channel+"_"+args.ftype+"_logy.png")

        c.Clear()
        pad.Delete()
        c.Delete()


        f_1= ss_1_tight.Clone()
        f_1.Divide(ss_1_loose)

        f_1.GetYaxis().SetTitleOffset(1.4)
        f_1.GetYaxis().SetMaxDigits(2)
        ROOT.TGaxis().SetMaxDigits(2)

        f_1.SetName("FakeRate_"+args.lepton)
        f_1.SetTitle(" Fake Rate Measurement "+args.lepton)
        f_1.GetXaxis().SetTitle("Leg p_T "+args.lepton)
        f_1.GetYaxis().SetTitle("Fake Rate for Leg"+args.lepton)

        tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax())

        c=ROOT.TCanvas("canvas","",0,0,600,600)
        ROOT.gStyle.SetOptStat(0)
        f_1.Draw()
        f_1.Fit("tf_1")
        c.SaveAs("FFhistos_"+str(args.ffout)+"/"+args.channel+"_"+args.lepton+".png")
        tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
        f_1.Write("fakerate1",ROOT.TObject.kOverwrite)

        ss_1_tight.Write("ss_1_tight_sub",ROOT.TObject.kOverwrite)
        ss_1_loose.Write("ss_1_loose_sub",ROOT.TObject.kOverwrite)




        datadrivenPackage["fakerate1"]=f_1.Clone()
        datadrivenPackage["fitrate1"]=tf_1
        datadrivenPackage["fakemeasurefile"]=fakemeasurefile
        fakemeasurefile.Write()


    if args.ftype=="muon":


        cattight_mt = "mmmt_FF_SS_validation"
        catloose_mt = "mmmt_FF_SS_1_loose"

        cattight_em = "mmem_FF_SS_validation"
        catloose_em = "mmem_FF_SS_2_loose"

        measurecats = [cattight_mt,cattight_em,catloose_mt,catloose_em]
        for cat in measurecats:
            infofile.write(cat)
            infofile.write("\n")
            #print(allcats[cat].cuts)
            for i,cut in allcats[cat].cuts.items():
                infofile.write(i)
                infofile.write(" : ")
                for subcut in cut:
                    infofile.write(str(subcut))
                infofile.write("\n")

        ss_1_tight_mt = histodict[cattight_mt]["Nominal_data_obs"]["pt_3_ff"]
        ss_1_loose_mt = histodict[catloose_mt]["Nominal_data_obs"]["pt_3_ff"]
        ss_1_tight_em = histodict[cattight_em]["Nominal_data_obs"]["pt_4_ff"]
        ss_1_loose_em = histodict[catloose_em]["Nominal_data_obs"]["pt_4_ff"]
        if args.prompt:
            ss_1_tight_mt_prompt = histodict[cattight_mt]["Nominal_prompt1"]["pt_3_ff"]
            ss_1_loose_mt_prompt = histodict[catloose_mt]["Nominal_prompt1"]["pt_3_ff"]
            ss_1_tight_em_prompt = histodict[cattight_em]["Nominal_prompt2"]["pt_4_ff"]
            ss_1_loose_em_prompt = histodict[catloose_em]["Nominal_prompt2"]["pt_4_ff"]

        c = ROOT.TCanvas()
        pad = ROOT.TPad()
        padlist=[]
        leg = ROOT.TLegend()
        doRatio=False
        c,padlist,leg = setUpPlot(doRatio)
        lumi,cms,pre = getPlotStuffs(args.year,args.channel)

        pad = padlist[0]
        ROOT.gStyle.SetOptStat(0)
        ###################
        #mutau
        ###################
        c.SetLogy(0)
        ss_1_tight_mt_fine=histodict[cattight_mt]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[cattight_mt]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Numerator")
        histodict[cattight_mt]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[cattight_mt]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[cattight_mt]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_mmmt_"+args.ftype+".png")
        c.SetLogy(1)
        histodict[cattight_mt]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_mmmt_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)

        ss_1_loose_mt_fine=histodict[catloose_mt]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[catloose_mt]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Denominator")
        histodict[catloose_mt]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[catloose_mt]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[catloose_mt]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_mmmt_"+args.ftype+".png")

        c.SetLogy(1)
        histodict[catloose_mt]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_mmmt_"+args.ftype+"_logy.png")

        c.Clear()
        ###################
        #ele muon
        ###################
        c.SetLogy(0)
        ss_1_tight_em_fine=histodict[cattight_em]["Nominal_data_obs"]["pt_4_ff_fine"]
        histodict[cattight_em]["Nominal_data_obs"]["pt_4_ff_fine"].SetTitle("Lepton p_T in Numerator")
        histodict[cattight_em]["Nominal_data_obs"]["pt_4_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[cattight_em]["Nominal_data_obs"]["pt_4_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[cattight_em]["Nominal_data_obs"]["pt_4_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_mmem_"+args.ftype+".png")
        c.SetLogy(1)
        histodict[cattight_em]["Nominal_data_obs"]["pt_4_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_mmem_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)

        ss_1_loose_em_fine=histodict[catloose_em]["Nominal_data_obs"]["pt_4_ff_fine"]
        histodict[catloose_em]["Nominal_data_obs"]["pt_4_ff_fine"].SetTitle("Lepton p_T in Denominator")
        histodict[catloose_em]["Nominal_data_obs"]["pt_4_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[catloose_em]["Nominal_data_obs"]["pt_4_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[catloose_em]["Nominal_data_obs"]["pt_4_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_mmem_"+args.ftype+".png")

        c.SetLogy(1)
        histodict[catloose_em]["Nominal_data_obs"]["pt_4_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_mmem_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)


        f_1_fine= ss_1_tight_mt_fine.Clone()
        f_1_fine.Add(ss_1_tight_em_fine)
        d_1_fine = ss_1_loose_mt_fine.Clone()
        d_1_fine.Add(ss_1_loose_em_fine)


        f_1_fine.SetTitle("Lepton p_T Numerator full #mu")
        f_1_fine.GetXaxis().SetTitle("Lepton p_T")
        f_1_fine.GetYaxis().SetTitle("Events ")
        f_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_full_"+args.ftype+".png")

        c.Clear()
        c.SetLogy(1)
        f_1_fine.SetTitle("Lepton p_T Numerator full #mu")
        f_1_fine.GetXaxis().SetTitle("Lepton p_T")
        f_1_fine.GetYaxis().SetTitle("Events ")
        f_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_full_"+args.ftype+"_logy.png")
        c.Clear()

        c.SetLogy(0)
        d_1_fine.SetTitle("Lepton p_T Denominator full #mu")
        d_1_fine.GetXaxis().SetTitle("Lepton p_T")
        d_1_fine.GetYaxis().SetTitle("Events ")
        d_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_full_"+args.ftype+".png")

        c.Clear()
        c.SetLogy(1)
        d_1_fine.SetTitle("Lepton p_T Denominator full #mu")
        d_1_fine.GetXaxis().SetTitle("Lepton p_T")
        d_1_fine.GetYaxis().SetTitle("Events ")
        d_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_full_"+args.ftype+"_logy.png")
        c.Clear()

        pad.Delete()
        c.Delete()

        f_1= ss_1_tight_mt.Clone()
        f_1.Add(ss_1_tight_em)
        d_1 = ss_1_loose_mt.Clone()
        d_1.Add(ss_1_loose_em)

        f_1.Divide(d_1)

        f_1.GetYaxis().SetTitleOffset(1.4)
        f_1.GetYaxis().SetMaxDigits(2)

        f_1.SetName("FakeRate_"+args.ftype)
        f_1.SetTitle(" Fake Rate Measurement "+args.ftype)
        f_1.GetXaxis().SetTitle("p_T "+args.ftype)
        f_1.GetYaxis().SetTitle("Fake Rate for "+args.ftype)

        tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax())

        c=ROOT.TCanvas("canvas","",0,0,600,600)
        ROOT.gStyle.SetOptStat(0)
        Ylow=0.335
        Xlow=0.58
        lumi=add_lumi(Xlow,Ylow,args.year,doRatio)
        Xlow=0.27
        cms=add_CMS(Xlow,Ylow,doRatio)
        Xlow=0.39
        pre=add_Preliminary(Xlow,Ylow,args.channel, doRatio)
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        f_1.Draw()
        f_1.Fit("tf_1")
        c.SaveAs("FFhistos_"+str(args.ffout)+"/"+args.channel+"_"+args.ftype+".png")
        tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
        f_1.Write("fakerate1",ROOT.TObject.kOverwrite)

        ss_1_tight_mt.Write("ss_1_tight_mt_sub",ROOT.TObject.kOverwrite)
        ss_1_loose_mt.Write("ss_1_loose_mt_sub",ROOT.TObject.kOverwrite)
        ss_1_tight_em.Write("ss_1_tight_em_sub",ROOT.TObject.kOverwrite)
        ss_1_loose_em.Write("ss_1_loose_em_sub",ROOT.TObject.kOverwrite)



        datadrivenPackage["fakerate1"]=f_1.Clone()
        datadrivenPackage["fitrate1"]=tf_1
        datadrivenPackage["fakemeasurefile"]=fakemeasurefile
        fakemeasurefile.Write()

        print("data in tight mutau ",ss_1_tight_mt.Integral())
        infofile.write("data in tight mutau\n")
        infofile.write(str(ss_1_tight_mt.Integral()))
        infofile.write("\n")
        print("data in tight electronmuon ",ss_1_tight_em.Integral())
        infofile.write("data in tight electronmuon  \n")
        infofile.write(str(ss_1_tight_em.Integral()))
        infofile.write("\n")

        print("data in loose mutau ",ss_1_loose_mt.Integral())
        infofile.write("data in loose mutau \n")
        infofile.write(str(ss_1_loose_mt.Integral()))
        infofile.write("\n")
        print("data in loose electronmuon ",ss_1_loose_em.Integral())
        infofile.write("data in loose electronmuon \n")
        infofile.write(str(ss_1_loose_em.Integral()))
        infofile.write("\n")

        if args.prompt:
            print("prompt in tight mutau ",ss_1_tight_mt_prompt.Integral())
            infofile.write("prompt in tight mutau\n")
            infofile.write(str(ss_1_tight_mt_prompt.Integral()))
            infofile.write("\n")
            print("prompt in tight electronmuon ",ss_1_tight_em_prompt.Integral())
            infofile.write("prompt in tight electronmuon  \n")
            infofile.write(str(ss_1_tight_em_prompt.Integral()))
            infofile.write("\n")

            print("prompt in loose mutau ",ss_1_loose_mt_prompt.Integral())
            infofile.write("prompt in loose mutau \n")
            infofile.write(str(ss_1_loose_mt_prompt.Integral()))
            infofile.write("\n")
            print("prompt in loose electronmuon ",ss_1_loose_em_prompt.Integral())
            infofile.write("prompt in loose electronmuon \n")
            infofile.write(str(ss_1_loose_em_prompt.Integral()))
            infofile.write("\n")

    if args.ftype=="electron":
        cattight_et = "mmet_FF_SS_validation"
        catloose_et = "mmet_FF_SS_1_loose"

        cattight_em = "mmem_FF_SS_validation"
        catloose_em = "mmem_FF_SS_1_loose"

        measurecats = [cattight_et,cattight_em,catloose_et,catloose_em]
        for cat in measurecats:
            infofile.write(cat)
            infofile.write("\n")
            #print(allcats[cat].cuts)
            for i,cut in allcats[cat].cuts.items():
                infofile.write(i)
                infofile.write(" : ")
                for subcut in cut:
                    infofile.write(str(subcut))
                infofile.write("\n")

        ss_1_tight_et = histodict[cattight_et]["Nominal_data_obs"]["pt_3_ff"]
        ss_1_loose_et = histodict[catloose_et]["Nominal_data_obs"]["pt_3_ff"]
        ss_1_tight_em = histodict[cattight_em]["Nominal_data_obs"]["pt_3_ff"]
        ss_1_loose_em = histodict[catloose_em]["Nominal_data_obs"]["pt_3_ff"]

        if args.prompt:
            ss_1_tight_et_prompt = histodict[cattight_et]["Nominal_prompt1"]["pt_3_ff"]
            ss_1_loose_et_prompt = histodict[catloose_et]["Nominal_prompt1"]["pt_3_ff"]
            ss_1_tight_em_prompt = histodict[cattight_em]["Nominal_prompt1"]["pt_3_ff"]
            ss_1_loose_em_prompt = histodict[catloose_em]["Nominal_prompt1"]["pt_3_ff"]

        c = ROOT.TCanvas()
        pad = ROOT.TPad()
        padlist=[]
        leg = ROOT.TLegend()
        doRatio=False
        c,padlist,leg = setUpPlot(doRatio)
        lumi,cms,pre = getPlotStuffs(args.year,args.channel)

        pad = padlist[0]
        ROOT.gStyle.SetOptStat(0)
        ###################
        #mutau
        ###################
        c.SetLogy(0)
        ss_1_tight_et_fine=histodict[cattight_et]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[cattight_et]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Numerator")
        histodict[cattight_et]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[cattight_et]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[cattight_et]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_mmet_"+args.ftype+".png")
        c.SetLogy(1)
        histodict[cattight_et]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_mmet_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)

        ss_1_loose_et_fine=histodict[catloose_et]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[catloose_et]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Denominator")
        histodict[catloose_et]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[catloose_et]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[catloose_et]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_mmet_"+args.ftype+".png")

        c.SetLogy(1)
        histodict[catloose_et]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_mmet_"+args.ftype+"_logy.png")

        c.Clear()
        ###################
        #ele muon
        ###################
        c.SetLogy(0)
        ss_1_tight_em_fine=histodict[cattight_em]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[cattight_em]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Numerator")
        histodict[cattight_em]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[cattight_em]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[cattight_em]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_mmem_"+args.ftype+".png")
        c.SetLogy(1)
        histodict[cattight_em]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_mmem_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)

        ss_1_loose_em_fine=histodict[catloose_em]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[catloose_em]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Denominator")
        histodict[catloose_em]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[catloose_em]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[catloose_em]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_mmem_"+args.ftype+".png")

        c.SetLogy(1)
        histodict[catloose_em]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_mmem_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)


        f_1_fine= ss_1_tight_et_fine.Clone()
        f_1_fine.Add(ss_1_tight_em_fine)
        d_1_fine = ss_1_loose_et_fine.Clone()
        d_1_fine.Add(ss_1_loose_em_fine)


        f_1_fine.SetTitle("Lepton p_T Numerator full e")
        f_1_fine.GetXaxis().SetTitle("Lepton p_T")
        f_1_fine.GetYaxis().SetTitle("Events ")
        f_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_full_"+args.ftype+".png")

        c.Clear()
        c.SetLogy(1)
        f_1_fine.SetTitle("Lepton p_T Numerator full e")
        f_1_fine.GetXaxis().SetTitle("Lepton p_T")
        f_1_fine.GetYaxis().SetTitle("Events ")
        f_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_full_"+args.ftype+"_logy.png")
        c.Clear()

        c.SetLogy(0)
        d_1_fine.SetTitle("Lepton p_T Denominator full e")
        d_1_fine.GetXaxis().SetTitle("Lepton p_T")
        d_1_fine.GetYaxis().SetTitle("Events ")
        d_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_full_"+args.ftype+".png")

        c.Clear()
        c.SetLogy(1)
        d_1_fine.SetTitle("Lepton p_T Denominator full e")
        d_1_fine.GetXaxis().SetTitle("Lepton p_T")
        d_1_fine.GetYaxis().SetTitle("Events ")
        d_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_full_"+args.ftype+"_logy.png")
        c.Clear()

        pad.Delete()
        c.Delete()

        f_1= ss_1_tight_et.Clone()
        f_1.Add(ss_1_tight_em)
        d_1 = ss_1_loose_et.Clone()
        d_1.Add(ss_1_loose_em)

        f_1.Divide(d_1)

        f_1.GetYaxis().SetTitleOffset(1.4)
        f_1.GetYaxis().SetMaxDigits(2)

        f_1.SetName("FakeRate_"+args.ftype)
        f_1.SetTitle(" Fake Rate Measurement "+args.ftype)
        f_1.GetXaxis().SetTitle("p_T "+args.ftype)
        f_1.GetYaxis().SetTitle("Fake Rate for "+args.ftype)

        tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax())

        c=ROOT.TCanvas("canvas","",0,0,600,600)
        ROOT.gStyle.SetOptStat(0)
        f_1.Draw()
        f_1.Fit("tf_1")
        Ylow=0.335
        Xlow=0.58
        lumi=add_lumi(Xlow,Ylow,args.year,doRatio)
        Xlow=0.27
        cms=add_CMS(Xlow,Ylow,doRatio)
        Xlow=0.39
        pre=add_Preliminary(Xlow,Ylow,args.channel, doRatio)
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/"+args.channel+"_"+args.ftype+".png")
        tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
        f_1.Write("fakerate1",ROOT.TObject.kOverwrite)

        ss_1_tight_et.Write("ss_1_tight_et_sub",ROOT.TObject.kOverwrite)
        ss_1_loose_et.Write("ss_1_loose_et_sub",ROOT.TObject.kOverwrite)
        ss_1_tight_em.Write("ss_1_tight_em_sub",ROOT.TObject.kOverwrite)
        ss_1_loose_em.Write("ss_1_loose_em_sub",ROOT.TObject.kOverwrite)



        datadrivenPackage["fakerate1"]=f_1.Clone()
        datadrivenPackage["fitrate1"]=tf_1
        datadrivenPackage["fakemeasurefile"]=fakemeasurefile
        fakemeasurefile.Write()

        print("data in tight electrontau ",ss_1_tight_et.Integral())
        infofile.write("data in tight electrontau\n")
        infofile.write(str(ss_1_tight_et.Integral()))
        infofile.write("\n")
        print("data in tight electronmuon ",ss_1_tight_em.Integral())
        infofile.write("data in tight electronmuon  \n")
        infofile.write(str(ss_1_tight_em.Integral()))
        infofile.write("\n")

        print("data in loose electrontau ",ss_1_loose_et.Integral())
        infofile.write("data in loose electrontau \n")
        infofile.write(str(ss_1_loose_et.Integral()))
        infofile.write("\n")
        print("data in loose electronmuon ",ss_1_loose_em.Integral())
        infofile.write("data in loose electronmuon \n")
        infofile.write(str(ss_1_loose_em.Integral()))
        infofile.write("\n")

        if args.prompt:
            print("prompt in tight electrontau ",ss_1_tight_et_prompt.Integral())
            infofile.write("prompt in tight electrontau\n")
            infofile.write(str(ss_1_tight_et_prompt.Integral()))
            infofile.write("\n")
            print("prompt in tight electronmuon ",ss_1_tight_em_prompt.Integral())
            infofile.write("prompt in tight electronmuon  \n")
            infofile.write(str(ss_1_tight_em_prompt.Integral()))
            infofile.write("\n")

            print("prompt in loose electrontau ",ss_1_loose_et_prompt.Integral())
            infofile.write("prompt in loose electrontau \n")
            infofile.write(str(ss_1_loose_et_prompt.Integral()))
            infofile.write("\n")
            print("prompt in loose electronmuon ",ss_1_loose_em_prompt.Integral())
            infofile.write("prompt in loose electronmuon \n")
            infofile.write(str(ss_1_loose_em_prompt.Integral()))
            infofile.write("\n")

    ##################
    # measuring tau with both legs
    ##################
    if args.ftype=="tautau":
        cattight_t = "mmtt_FF_SS_validation"
        catloose_1 = "mmtt_FF_SS_1_loose"
        catloose_2 = "mmtt_FF_SS_2_loose"

        measurecats = [cattight_t,catloose_1,catloose_2]
        for cat in measurecats:
            infofile.write(cat)
            infofile.write("\n")
            #print(allcats[cat].cuts)
            for i,cut in allcats[cat].cuts.items():
                infofile.write(i)
                infofile.write(" : ")
                for subcut in cut:
                    infofile.write(str(subcut))
                infofile.write("\n")

        ss_1_tight_t = histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff"]
        ss_1_loose_1 = histodict[catloose_1]["Nominal_data_obs"]["pt_3_ff"]
        ss_1_loose_2 = histodict[catloose_2]["Nominal_data_obs"]["pt_4_ff"]
        if args.prompt:
            ss_1_tight_t_prompt = histodict[cattight_t]["Nominal_prompt1"]["pt_3_ff"]
            ss_1_loose_1_prompt = histodict[catloose_1]["Nominal_prompt1"]["pt_3_ff"]
            ss_1_loose_2_prompt = histodict[catloose_2]["Nominal_prompt2"]["pt_4_ff"]

        c = ROOT.TCanvas()
        pad = ROOT.TPad()
        padlist=[]
        leg = ROOT.TLegend()
        doRatio=False
        c,padlist,leg = setUpPlot(doRatio)
        lumi,cms,pre = getPlotStuffs(args.year,args.channel)

        pad = padlist[0]
        ROOT.gStyle.SetOptStat(0)
        ###################
        #tau
        ###################
        c.SetLogy(0)
        ss_1_tight_t_fine=histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Numerator")
        histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_tau_"+args.ftype+".png")
        c.SetLogy(1)
        histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_tau_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)

        ss_1_loose_1_fine=histodict[catloose_1]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[catloose_1]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Denominator")
        histodict[catloose_1]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[catloose_1]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[catloose_1]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_tau_"+args.ftype+".png")

        c.SetLogy(1)
        histodict[catloose_1]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_tau_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)

        ss_1_loose_2_fine=histodict[catloose_2]["Nominal_data_obs"]["pt_4_ff_fine"]
        histodict[catloose_2]["Nominal_data_obs"]["pt_4_ff_fine"].SetTitle("Lepton p_T in Denominator")
        histodict[catloose_2]["Nominal_data_obs"]["pt_4_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[catloose_2]["Nominal_data_obs"]["pt_4_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[catloose_2]["Nominal_data_obs"]["pt_4_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_muon_"+args.ftype+".png")

        c.SetLogy(1)
        histodict[catloose_2]["Nominal_data_obs"]["pt_4_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_muon_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)


        f_1_fine= ss_1_tight_t_fine.Clone()
        f_1_fine.Add(ss_1_tight_t_fine)
        d_1_fine = ss_1_loose_1_fine.Clone()
        d_1_fine.Add(ss_1_loose_2_fine)


        f_1_fine.SetTitle("Lepton p_T Numerator full #tau")
        f_1_fine.GetXaxis().SetTitle("Lepton p_T")
        f_1_fine.GetYaxis().SetTitle("Events ")
        f_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_full_"+args.ftype+".png")

        c.Clear()
        c.SetLogy(1)
        f_1_fine.SetTitle("Lepton p_T Numerator full #tau")
        f_1_fine.GetXaxis().SetTitle("Lepton p_T")
        f_1_fine.GetYaxis().SetTitle("Events ")
        f_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_full_"+args.ftype+"_logy.png")
        c.Clear()

        c.SetLogy(0)
        d_1_fine.SetTitle("Lepton p_T Denominator full #tau")
        d_1_fine.GetXaxis().SetTitle("Lepton p_T")
        d_1_fine.GetYaxis().SetTitle("Events ")
        d_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_full_"+args.ftype+".png")

        c.Clear()
        c.SetLogy(1)
        d_1_fine.SetTitle("Lepton p_T Denominator full #tau")
        d_1_fine.GetXaxis().SetTitle("Lepton p_T")
        d_1_fine.GetYaxis().SetTitle("Events ")
        d_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_full_"+args.ftype+"_logy.png")
        c.Clear()

        pad.Delete()
        c.Delete()

        f_1= ss_1_tight_t.Clone()
        f_1.Add(ss_1_tight_t)
        d_1 = ss_1_loose_1.Clone()
        d_1.Add(ss_1_loose_2)

        f_1.Divide(d_1)

        f_1.GetYaxis().SetTitleOffset(1.4)
        f_1.GetYaxis().SetMaxDigits(2)

        f_1.SetName("FakeRate"+args.ftype)
        f_1.SetTitle(" Fake Rate Measurement "+args.ftype)
        f_1.GetXaxis().SetTitle("p_T "+args.ftype)
        f_1.GetYaxis().SetTitle("Fake Rate for "+args.ftype)

        tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax())

        c=ROOT.TCanvas("canvas","",0,0,600,600)
        ROOT.gStyle.SetOptStat(0)
        f_1.Draw()
        f_1.Fit("tf_1")
        Ylow=0.335
        Xlow=0.58
        lumi=add_lumi(Xlow,Ylow,args.year,doRatio)
        Xlow=0.27
        cms=add_CMS(Xlow,Ylow,doRatio)
        Xlow=0.39
        pre=add_Preliminary(Xlow,Ylow,args.channel, doRatio)
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/"+args.channel+"_"+args.ftype+".png")
        tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
        f_1.Write("fakerate1",ROOT.TObject.kOverwrite)

        ss_1_tight_t.Write("ss_1_tight_t_sub",ROOT.TObject.kOverwrite)
        ss_1_loose_1.Write("ss_1_loose_1_sub",ROOT.TObject.kOverwrite)
        ss_1_loose_2.Write("ss_1_loose_2_sub",ROOT.TObject.kOverwrite)



        datadrivenPackage["fakerate1"]=f_1.Clone()
        datadrivenPackage["fitrate1"]=tf_1
        datadrivenPackage["fakemeasurefile"]=fakemeasurefile
        fakemeasurefile.Write()

        print("data in tight tau ",ss_1_tight_t.Integral())
        infofile.write("data in tight tau\n")
        infofile.write(str(ss_1_tight_t.Integral()))
        infofile.write("\n")

        print("data in loose tau ",ss_1_loose_1.Integral())
        infofile.write("data in loose tau \n")
        infofile.write(str(ss_1_loose_1.Integral()))
        infofile.write("\n")
        print("data in loose muon ",ss_1_loose_2.Integral())
        infofile.write("data in loose muon \n")
        infofile.write(str(ss_1_loose_2.Integral()))
        infofile.write("\n")
        if args.prompt:
            print("prompt in tight tau ",ss_1_tight_t_prompt.Integral())
            infofile.write("prompt in tight tau\n")
            infofile.write(str(ss_1_tight_t_prompt.Integral()))
            infofile.write("\n")

            print("prompt in loose tau ",ss_1_loose_1_prompt.Integral())
            infofile.write("prompt in loose tau \n")
            infofile.write(str(ss_1_loose_1_prompt.Integral()))
            infofile.write("\n")
            print("prompt in loose tau 2 ",ss_1_loose_2_prompt.Integral())
            infofile.write("prompt in loose tau 2  \n")
            infofile.write(str(ss_1_loose_2_prompt.Integral()))
            infofile.write("\n")

    ##################
    # measuring tau with all channels  
    ##################
    if args.ftype=="tau":
        cattight_t = "mmtt_FF_SS_validation"
        cattight_m = "mmmt_FF_SS_validation"
        cattight_e = "mmet_FF_SS_validation"

        catloose_t = "mmtt_FF_SS_1_loose"
        catloose_m = "mmmt_FF_SS_2_loose"
        catloose_e = "mmet_FF_SS_2_loose"

        measurecats = [cattight_t,cattight_m,cattight_e,catloose_t,catloose_m,catloose_e]
        for cat in measurecats:
            infofile.write(cat)
            infofile.write("\n")
            #print(allcats[cat].cuts)
            for i,cut in allcats[cat].cuts.items():
                infofile.write(i)
                infofile.write(" : ")
                for subcut in cut:
                    infofile.write(str(subcut))
                infofile.write("\n")

        ss_1_tight_t = histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff"]
        ss_1_loose_t = histodict[catloose_t]["Nominal_data_obs"]["pt_3_ff"]
        ss_1_tight_m = histodict[cattight_m]["Nominal_data_obs"]["pt_4_ff"]
        ss_1_loose_m = histodict[catloose_m]["Nominal_data_obs"]["pt_4_ff"]
        ss_1_tight_e = histodict[cattight_e]["Nominal_data_obs"]["pt_4_ff"]
        ss_1_loose_e = histodict[catloose_e]["Nominal_data_obs"]["pt_4_ff"]
        if args.prompt:
            ss_1_tight_t_prompt = histodict[cattight_t]["Nominal_prompt1"]["pt_3_ff"]
            ss_1_loose_t_prompt = histodict[catloose_t]["Nominal_prompt1"]["pt_3_ff"]
            ss_1_tight_m_prompt = histodict[cattight_m]["Nominal_prompt2"]["pt_4_ff"]
            ss_1_loose_m_prompt = histodict[catloose_m]["Nominal_prompt2"]["pt_4_ff"]
            ss_1_tight_e_prompt = histodict[cattight_e]["Nominal_prompt2"]["pt_4_ff"]
            ss_1_loose_e_prompt = histodict[catloose_e]["Nominal_prompt2"]["pt_4_ff"]

        c = ROOT.TCanvas()
        pad = ROOT.TPad()
        padlist=[]
        leg = ROOT.TLegend()
        doRatio=False
        c,padlist,leg = setUpPlot(doRatio)
        lumi,cms,pre = getPlotStuffs(args.year,args.channel)

        pad = padlist[0]
        ROOT.gStyle.SetOptStat(0)
        ###################
        #tau
        ###################
        c.SetLogy(0)
        ss_1_tight_t_fine=histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Numerator")
        histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_tau_"+args.ftype+".png")
        c.SetLogy(1)
        histodict[cattight_t]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_tau_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)

        ss_1_loose_t_fine=histodict[catloose_t]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[catloose_t]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Denominator")
        histodict[catloose_t]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[catloose_t]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[catloose_t]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_tau_"+args.ftype+".png")

        c.SetLogy(1)
        histodict[catloose_t]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_tau_"+args.ftype+"_logy.png")

        c.Clear()
        ###################
        #muon
        ###################
        c.SetLogy(0)
        ss_1_tight_m_fine=histodict[cattight_m]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[cattight_m]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Numerator")
        histodict[cattight_m]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[cattight_m]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[cattight_m]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_muon_"+args.ftype+".png")
        c.SetLogy(1)
        histodict[cattight_m]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_muon_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)

        ss_1_loose_m_fine=histodict[catloose_m]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[catloose_m]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Denominator")
        histodict[catloose_m]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[catloose_m]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[catloose_m]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_muon_"+args.ftype+".png")

        c.SetLogy(1)
        histodict[catloose_m]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_muon_"+args.ftype+"_logy.png")

        ###################
        #electron
        ###################
        c.SetLogy(0)
        ss_1_tight_e_fine=histodict[cattight_e]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[cattight_e]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Numerator")
        histodict[cattight_e]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[cattight_e]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[cattight_e]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_electron_"+args.ftype+".png")
        c.SetLogy(1)
        histodict[cattight_e]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_electron_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)

        ss_1_loose_e_fine=histodict[catloose_e]["Nominal_data_obs"]["pt_3_ff_fine"]
        histodict[catloose_e]["Nominal_data_obs"]["pt_3_ff_fine"].SetTitle("Lepton p_T in Denominator")
        histodict[catloose_e]["Nominal_data_obs"]["pt_3_ff_fine"].GetXaxis().SetTitle("Lepton p_T ")
        histodict[catloose_e]["Nominal_data_obs"]["pt_3_ff_fine"].GetYaxis().SetTitle("Events ")
        histodict[catloose_e]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_electron_"+args.ftype+".png")

        c.SetLogy(1)
        histodict[catloose_e]["Nominal_data_obs"]["pt_3_ff_fine"].Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_electron_"+args.ftype+"_logy.png")

        c.Clear()

        c.SetLogy(0)


        f_1_fine= ss_1_tight_t_fine.Clone()
        f_1_fine.Add(ss_1_tight_m_fine)
        f_1_fine.Add(ss_1_tight_e_fine)
        d_1_fine = ss_1_loose_t_fine.Clone()
        d_1_fine.Add(ss_1_loose_m_fine)
        d_1_fine.Add(ss_1_loose_e_fine)


        f_1_fine.SetTitle("Lepton p_T Numerator full #tau")
        f_1_fine.GetXaxis().SetTitle("Lepton p_T")
        f_1_fine.GetYaxis().SetTitle("Events ")
        f_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_full_"+args.ftype+".png")

        c.Clear()
        c.SetLogy(1)
        f_1_fine.SetTitle("Lepton p_T Numerator full #tau")
        f_1_fine.GetXaxis().SetTitle("Lepton p_T")
        f_1_fine.GetYaxis().SetTitle("Events ")
        f_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Numerator_"+args.channel+"_full_"+args.ftype+"_logy.png")
        c.Clear()

        c.SetLogy(0)
        d_1_fine.SetTitle("Lepton p_T Denominator full #tau")
        d_1_fine.GetXaxis().SetTitle("Lepton p_T")
        d_1_fine.GetYaxis().SetTitle("Events ")
        d_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_full_"+args.ftype+".png")

        c.Clear()
        c.SetLogy(1)
        d_1_fine.SetTitle("Lepton p_T Denominator full #tau")
        d_1_fine.GetXaxis().SetTitle("Lepton p_T")
        d_1_fine.GetYaxis().SetTitle("Events ")
        d_1_fine.Draw("PE")
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/leptonPt_Denominator_"+args.channel+"_full_"+args.ftype+"_logy.png")
        c.Clear()

        pad.Delete()
        c.Delete()

        f_1= ss_1_tight_t.Clone()
        f_1.Add(ss_1_tight_m)
        f_1.Add(ss_1_tight_e)
        d_1 = ss_1_loose_t.Clone()
        d_1.Add(ss_1_loose_m)
        d_1.Add(ss_1_loose_e)

        f_1.Divide(d_1)

        f_1.GetYaxis().SetTitleOffset(1.4)
        f_1.GetYaxis().SetMaxDigits(2)

        f_1.SetName("FakeRate"+args.ftype)
        f_1.SetTitle(" Fake Rate Measurement "+args.ftype)
        f_1.GetXaxis().SetTitle("p_T "+args.ftype)
        f_1.GetYaxis().SetTitle("Fake Rate for "+args.ftype)

        tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax())

        c=ROOT.TCanvas("canvas","",0,0,600,600)
        ROOT.gStyle.SetOptStat(0)
        f_1.Draw()
        f_1.Fit("tf_1")
        Ylow=0.335
        Xlow=0.58
        lumi=add_lumi(Xlow,Ylow,args.year,doRatio)
        Xlow=0.27
        cms=add_CMS(Xlow,Ylow,doRatio)
        Xlow=0.39
        pre=add_Preliminary(Xlow,Ylow,args.channel, doRatio)
        lumi.Draw()
        cms.Draw()
        pre.Draw()
        c.SaveAs("FFhistos_"+str(args.ffout)+"/"+args.channel+"_"+args.ftype+".png")
        tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
        f_1.Write("fakerate1",ROOT.TObject.kOverwrite)

        ss_1_tight_t.Write("ss_1_tight_t_sub",ROOT.TObject.kOverwrite)
        ss_1_loose_t.Write("ss_1_loose_t_sub",ROOT.TObject.kOverwrite)
        ss_1_tight_m.Write("ss_1_tight_m_sub",ROOT.TObject.kOverwrite)
        ss_1_loose_m.Write("ss_1_loose_m_sub",ROOT.TObject.kOverwrite)
        ss_1_tight_e.Write("ss_1_tight_e_sub",ROOT.TObject.kOverwrite)
        ss_1_loose_e.Write("ss_1_loose_e_sub",ROOT.TObject.kOverwrite)



        datadrivenPackage["fakerate1"]=f_1.Clone()
        datadrivenPackage["fitrate1"]=tf_1
        datadrivenPackage["fakemeasurefile"]=fakemeasurefile
        fakemeasurefile.Write()

        print("data in tight tau ",ss_1_tight_t.Integral())
        infofile.write("data in tight tau\n")
        infofile.write(str(ss_1_tight_t.Integral()))
        infofile.write("\n")
        print("data in tight muon ",ss_1_tight_m.Integral())
        infofile.write("data in tight muon  \n")
        infofile.write(str(ss_1_tight_m.Integral()))
        infofile.write("\n")
        print("data in tight electron ",ss_1_tight_e.Integral())
        infofile.write("data in tight electron \n")
        infofile.write(str(ss_1_tight_e.Integral()))
        infofile.write("\n")

        print("data in loose tau ",ss_1_loose_t.Integral())
        infofile.write("data in loose tau \n")
        infofile.write(str(ss_1_loose_t.Integral()))
        infofile.write("\n")
        print("data in loose muon ",ss_1_loose_m.Integral())
        infofile.write("data in loose muon \n")
        infofile.write(str(ss_1_loose_m.Integral()))
        infofile.write("\n")
        print("data in loose electron ",ss_1_loose_e.Integral())
        infofile.write("data in loose electron \n")
        infofile.write(str(ss_1_loose_e.Integral()))
        infofile.write("\n")
        if args.prompt:
            print("prompt in tight tau ",ss_1_tight_t_prompt.Integral())
            infofile.write("prompt in tight tau\n")
            infofile.write(str(ss_1_tight_t_prompt.Integral()))
            infofile.write("\n")
            print("prompt in tight muon ",ss_1_tight_m_prompt.Integral())
            infofile.write("prompt in tight muon  \n")
            infofile.write(str(ss_1_tight_m_prompt.Integral()))
            infofile.write("\n")
            print("prompt in tight electron ",ss_1_tight_e_prompt.Integral())
            infofile.write("prompt in tight electron \n")
            infofile.write(str(ss_1_tight_e_prompt.Integral()))
            infofile.write("\n")

            print("prompt in loose tau ",ss_1_loose_t_prompt.Integral())
            infofile.write("prompt in loose tau \n")
            infofile.write(str(ss_1_loose_t_prompt.Integral()))
            infofile.write("\n")
            print("prompt in loose muon ",ss_1_loose_m_prompt.Integral())
            infofile.write("prompt in loose muon \n")
            infofile.write(str(ss_1_loose_m_prompt.Integral()))
            infofile.write("\n")
            print("prompt in loose electron ",ss_1_loose_e_prompt.Integral())
            infofile.write("prompt in loose electron \n")
            infofile.write(str(ss_1_loose_e_prompt.Integral()))
            infofile.write("\n")

    try:
        datadrivenPackage["fakemeasurefile"].Close()
    except:
        print("problem closing fake measurement helper file ... shutting down anyway")

    print("computation time")
    print(datetime.datetime.now() - begin_time)
    print("arguments used")
    print(args)
