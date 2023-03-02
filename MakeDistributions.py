# !/usr/bin/env python
#########################
#Author: Sam Higginbotham
'''

* File Name : MakeDistributions.py

* Purpose : For visualizing datasets. Will create ultra skimmed root files or histograms for visualization and for limit setting later

* Usage:
python MakeDistributions.py -c cat_mmtt_2016.yaml  -csv MCsamples_2016_v6_yaml.csv  -i /path/to/input/root/files/ -p processes_special_mmtt.yaml -o main_outputdirectory -fo fake_factor_measure_outputdirectory -ch mmtt -s

* Creation Date : june-20-2020

* Last Modified : feb-16-2021

'''
#########################
import sys
import os
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kFatal
import uproot
import pandas
import json
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
import argparse
parser = argparse.ArgumentParser(description="This file generates root files containing Histograms ... files in utils contain selections and settings")
parser.add_argument("-o",  "--outname", default="",  help="postfix string")
parser.add_argument("-fi",  "--ffin", default="",  help="fake factor files")
parser.add_argument("-fi1",  "--ffin1", default="",  help="fake factor split to leg 1")
parser.add_argument("-fi2",  "--ffin2", default="",  help="fake factor split to leg 2")
parser.add_argument("-catt",  "--cattight", default="",  help="if measuring single leg pick one tight id/iso category")
parser.add_argument("-catl",  "--catloose", default="",  help="if measuring single leg pick one loose id/iso category")
parser.add_argument("-ft",  "--ftype", default="",  help="")
parser.add_argument("-lt",  "--lepton", default="",  help="")
parser.add_argument("-year",  "--year", default="2016",  help="Year")
parser.add_argument("-fo",  "--ffout", default="",  help="fake factor files to output")
parser.add_argument("-c",  "--categories", default="categories_array.yaml",  help="categories yaml file")
parser.add_argument("-ch",  "--channel", default="mmmt",  help="Please list the channel for fake factor histograms")
parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
parser.add_argument("-i",  "--dir", default="/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov7_basic_10_6_4/src/2016_v7/",  help="Input files")
parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
parser.add_argument("-dm",  "--datameasure", default=False,action='store_true',  help="Use DataDriven Method measure part")
parser.add_argument("-ps",  "--promptsub", default=False,action='store_true',  help="When making fake rate histograms, subtract the prompt contribution")
parser.add_argument("-dmOD",  "--onlydata", default=False,action='store_true',  help="only produced distributions derived from data ntuple.")
parser.add_argument("-dbg",  "--debug", default=False,action='store_true',  help="Disable the parallel processing mode")
parser.add_argument("-rp",  "--runprocesses", default="",  help="run only on selected processes input file csv format or txt ... delimeter is new line")
parser.add_argument("-sys",  "--systematics", default=False,action='store_true',  help="Run on systematic trees, don't run nominal at the same time")
parser.add_argument("-ddHAA",  "--datadrivenHAA", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-ddSM",  "--datadrivenSM", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-ex",  "--extract", default=False,action='store_true',  help="Additional Cuts for Extraction")
parser.add_argument("-ff",  "--makeFakeHistos", default=False,action='store_true',  help="Just make fake rate histos")
parser.add_argument("-v",  "--verbose", default=False,action='store_true',  help="print per event")
parser.add_argument("-t",  "--test", default=False,action='store_true',  help="only do 1 event to test code")
parser.add_argument("-s",  "--skim", default=False,action='store_true',  help="skim input files to make more TTrees")
parser.add_argument("-mt",  "--mt", default=False,action='store_true',  help="Use Multithreading")
parser.add_argument("-pt",  "--maxprint", default=False,action='store_true',  help="Print Info on cats and processes")
parser.add_argument("-combine",  "--combine", default=False,action='store_true',  help="just combine root files")
parser.add_argument("-rr",  "--rerun", default=False,action='store_true',  help="check mass output and run on remaining files")
parser.add_argument("-pc",  "--processingcores", default=12,type=int,help="Number of cores for multiprocessing")
args = parser.parse_args()


year = args.year

basedir=os.getenv("CMSSW_BASE")
sys.path.append(basedir)

#Setting up operators for cut string iterator
ops = { "==": operator.eq, "!=": operator.ne, ">": operator.gt, "<": operator.lt, ">=": operator.ge, "<=": operator.le, "band": operator.and_,"bor":operator.or_}


def subArrayLogic(evt,array):
    boo=ops[array[1]](returnValue(evt,array[0]),array[2])
    return boo


def returnArray(masterArray,variable):
   #if variable in ["njets","jpt_1","jeta_1","jpt_2","jeta_2","bpt_1","bpt_2","nbtag","beta_1","beta_2"]:
   #    val = masterArray[variable][:,0]
   #else:
   #print "working on branch ",variable
   #if len(masterArray[variable].shape) == 1:   # flat array important for jagged arrays of input data
   #    val = masterArray[variable]
   if "[" in variable:
       basevar = variable.split("[")[0]
       index = int(variable.split("[")[1].split("]")[0])
       val = masterArray[basevar][:,index]
   else:
       val = masterArray[variable]
   return val

def makeCutstar(package):
    return makeCut(*package)

def makeCut(mask,masterArray,cut):
    #print("making cut ",cut)
    if cut[0][0]=="SPECIAL":
        mask *= cutOnArray(masterArray,cut[1])
        return mask

    if cut[0][0]=="EQT":
        for i,var in enumerate(cut[1]):
            if cut[2]=="mult":
                if(i==0):
                    #tempmask=np.full(len(masterArray["evt"]),1.0)
                    tempmask=np.full(len(masterArray["mll"]),1.0)
                tempmask = tempmask * returnArray(masterArray,var)

            if cut[2]=="div":
                if(i==0):
                    #tempmask=np.full(len(masterArray["evt"]),1.0)
                    tempmask=np.full(len(masterArray["mll"]),1.0)
                tempmask = tempmask / returnArray(masterArray,var)
            if cut[2]=="add":
                if(i==0):
                    #tempmask=np.full(len(masterArray["evt"]),0.0)
                    tempmask=np.full(len(masterArray["mll"]),0.0)
                tempmask = tempmask + returnArray(masterArray,var)
            if cut[2]=="sub":
                if(i==0):
                    #tempmask=np.full(len(masterArray["evt"]),0.0)
                    tempmask=np.full(len(masterArray["mll"]),0.0)
                tempmask = tempmask - returnArray(masterArray,var)

        tempmask = ops[cut[3]](tempmask,cut[4])
        #test

        mask *= tempmask.astype(bool)
        return mask

    if cut[0][0]=="OR":
        tempmasks=[]
        #tempmask=np.full(len(masterArray["evt"]),False)
        tempmask=np.full(len(masterArray["mll"]),False)
        for oe in range(1,len(cut)):
            oneor = ops[cut[oe][1]](returnArray(masterArray,cut[oe][0]),cut[oe][2])
            #tempmask = ops["bor"](tempmask,oneor).astype(bool)
            tempmask = (tempmask + oneor).astype(bool)
        mask *= tempmask.astype(bool)
        return mask

    else:
        #print(cut)
        statement = cut[1]
        if statement=="absg": # lessthan OR greater than
            tempmask = ops["<"](returnArray(masterArray,cut[0]),(-1 * cut[2]))
            tempmask2= ops[">"](returnArray(masterArray,cut[0]),cut[2])
            tempmask = ops["bor"](tempmask,tempmask2)
        if statement=="absl": # greater than AND less tahn
            tempmask = ops[">"](returnArray(masterArray,cut[0]),(-1 * cut[2]))
            tempmask *= ops["<"](returnArray(masterArray,cut[0]),cut[2])
        else:
            tempmask = ops[statement](returnArray(masterArray,cut[0]),cut[2])
        mask *= tempmask.astype(bool)
        return mask
    return(mask)

def cutOnArray(masterArray,cuts):

    #mask=np.full(len(masterArray["evt"]),True)
    mask=np.full(len(masterArray["mll"]),True)

    packages = []
    masks = []
    for cut in cuts:
        packages.append((mask,masterArray,cut))
    #i = 0
    #m = mp.Manager()
    #logger_q = m.Queue()
    #parallelable_data = [(1, logger_q), (2, logger_q)]
    #pool  = mp.Pool(12)
    ##masks = pool.map(makeCutstar,packages)
    #pool.close()
    #pool.join()
        
    for pk in packages:
        masks.append(makeCutstar(pk))

    for mk in masks:
        mask*=mk

    return mask

from correctionlib import _core #requires cmssw 12_3_1  and python3 

basedir=os.getenv("CMSSW_BASE")
jsonpogloc=os.environ["CMSSW_BASE"]+"/src/jsonpog-integration/"
if year=="2016pre":
    tauFile = jsonpogloc+"POG/TAU/2016preVFP_UL/tau.json.gz"
    muFile = jsonpogloc+"POG/MUO/2016preVFP_UL/muon_Z.json.gz"
    eleFile = jsonpogloc+"POG/EGM/2016preVFP_UL/electron.json.gz"
if year=="2016post":
    tauFile = jsonpogloc+"POG/TAU/2016postVFP_UL/tau.json.gz"
    muFile = jsonpogloc+"POG/MUO/2016postVFP_UL/muon_Z.json.gz"
    eleFile = jsonpogloc+"POG/EGM/2016postVFP_UL/electron.json.gz"
if year=="2017":
    tauFile = jsonpogloc+"POG/TAU/2017_UL/tau.json.gz"
    muFile = jsonpogloc+"POG/MUO/2017_UL/muon_Z.json.gz"
    eleFile = jsonpogloc+"POG/EGM/2017_UL/electron.json.gz"
if year=="2018":
    tauFile = jsonpogloc+"POG/TAU/2018_UL/tau.json.gz"
    muFile = jsonpogloc+"POG/MUO/2018_UL/muon_Z.json.gz"
    eleFile = jsonpogloc+"POG/EGM/2018_UL/electron.json.gz"
cset_tau = _core.CorrectionSet.from_file(tauFile)
antiJetSF = cset_tau['DeepTau2017v2p1VSjet'] #'Loose'
antiEleSF = cset_tau['DeepTau2017v2p1VSe']  #'VLoose'
antiMuSF  = cset_tau['DeepTau2017v2p1VSmu'] #'Tight'

cset_ele = _core.CorrectionSet.from_file(eleFile)

def getAntiJetSF(pt,dm):
    #print(pt,dm)
    if args.channel == "mmmt":
        sf = antiJetSF.evaluate(float(pt),dm,5,"Loose","nom","pt")
    if args.channel == "mmet":
        sf = antiJetSF.evaluate(float(pt),dm,5,"Medium","nom","pt")
    if args.channel == "mmtt":
        sf = antiJetSF.evaluate(float(pt),dm,5,"VLoose","nom","pt")
    if args.channel == "mmem":
        sf = antiJetSF.evaluate(float(pt),dm,5,"Loose","nom","pt")
    if sf==0:
        sf=1.0
    return sf 

def getAntiEleSF(eta,gen_match):
    re_gm = 0
    if gen_match==15: re_gm=1
    if gen_match== 5: re_gm=5
    if gen_match==-1: re_gm=0
    if args.channel == "mmmt":
        sf = antiEleSF.evaluate(float(eta),re_gm,"VLoose","nom")
    if args.channel == "mmet":
        sf = antiEleSF.evaluate(float(eta),re_gm,"Medium","nom")
    if args.channel == "mmtt":
        sf = antiEleSF.evaluate(float(eta),re_gm,"VLoose","nom")
    if args.channel == "mmem":
        sf = antiEleSF.evaluate(float(eta),re_gm,"Medium","nom")
    if sf==0:
        sf=1.0
    return sf 

def getAntiMuSF(eta,gen_match):
    re_gm = 0
    if gen_match==15: re_gm=2
    if gen_match== 5: re_gm=5
    if gen_match==-1: re_gm=0
    #genmatch: 0 or 6 = unmatched or jet, 1 or 3 = electron, 2 or 4 = muon, 5 = real tau
    if args.channel == "mmmt":
        sf = antiMuSF.evaluate(float(eta),re_gm,"Tight","nom")
    if args.channel == "mmet":
        sf = antiMuSF.evaluate(float(eta),re_gm,"VLoose","nom")
    if args.channel == "mmtt":
        sf = antiMuSF.evaluate(float(eta),re_gm,"VLoose","nom")
    if args.channel == "mmem":
        sf = antiMuSF.evaluate(float(eta),re_gm,"Tight","nom")
    if sf==0:
        sf=1.0
    return sf 

#electronSFfile = ROOT.TFile()
year = args.year 
#if year=="2016":
#    electronSFfile = ROOT.TFile.Open(os.environ["CMSSW_BASE"]+"/src/AnalysisVisualization/2016LegacyReReco_ElectronMVA90noiso_Fall17V2.root")
#if year=="2017":
#    electronSFfile = ROOT.TFile.Open(os.environ["CMSSW_BASE"]+"/src/AnalysisVisualization/2017_ElectronMVA90noiso.root")
#if year=="2018":
#    electronSFfile = ROOT.TFile.Open(os.environ["CMSSW_BASE"]+"/src/AnalysisVisualization/2018_ElectronMVA90noiso.root")

#fact = self.weights_eleES.evaluate(year,esys,"wp90noiso",entry.Electron_eta[j],entry.Electron_pt[j])
eleES = cset_ele["UL-Electron-ID-SF"]
def getElectronSF(pt,eta):
    eyear = year
    if pt<10.0: return 1.0
    if "post" in year: eyear="2016postVFP"
    if "pre" in year: eyear="2016preVFP"
    sf = eleES.evaluate(eyear,"sf","wp90noiso",float(eta),float(pt))
    if sf==0:
        sf=1.0
    #print("Electron pt ",pt," Electron sf ",sf)
    return sf 

lowmuSFfile = ROOT.TFile()
muSFhisto = ROOT.TH2D()
#if year=="2016":
#    lowmuSFfile = ROOT.TFile.Open(os.environ["CMSSW_BASE"]+"/src/muonefficiencies/Run2/preUL/2016/2016_Jpsi/RunBCDEF_SF_ID.root")
#if year=="2017":
#    lowmuSFfile = ROOT.TFile.Open(os.environ["CMSSW_BASE"]+"/src/muonefficiencies/Run2/preUL/2017/2017_Jpsi/RunBCDEF_SF_ID.root")
#if year=="2018":
#    lowmuSFfile = ROOT.TFile.Open(os.environ["CMSSW_BASE"]+"/src/muonefficiencies/Run2/preUL/2018/2018_Jpsi/RunBCDEF_SF_ID.root")
if year=="2016pre":
    f = open(os.environ["CMSSW_BASE"]+"/src/muonefficiencies/Run2/UL/2016_preVFP/2016_preVFP_Jpsi/Efficiency_muon_trackerMuon_Run2016preVFP_UL_ID.json")
    sfs = json.load(f)

if year=="2016post":
    f = open(os.environ["CMSSW_BASE"]+"/src/muonefficiencies/Run2/UL/2016_postVFP/2016_postVFP_Jpsi/Efficiency_muon_trackerMuon_Run2016postVFP_UL_ID.json")
    sfs = json.load(f)

if year=="2017":

    f = open(os.environ["CMSSW_BASE"]+"/src/muonefficiencies/Run2/UL/2017/2017_Jpsi/Efficiency_muon_trackerMuon_Run2017_UL_ID.json")
    sfs = json.load(f)

if year=="2018":
    f = open(os.environ["CMSSW_BASE"]+"/src/muonefficiencies/Run2/UL/2018/2018_Jpsi/Efficiency_muon_trackerMuon_Run2018_UL_ID.json")
    sfs = json.load(f)

def getLowMuSF(pt,eta):
    sel_eta = 0.0 
    sel_pt = 0.0 
    sf = 0.0
    if pt<5.0: return 1.0
    for etarange in list(sfs["NUM_MediumID_DEN_TrackerMuons"]["abseta_pt"].keys()):
        if eta >= float(etarange.split("[")[1].split(",")[0]) and eta <= float(etarange.split("]")[0].split(",")[1]):
            sel_eta = etarange

    if sel_eta == 0.0: return 1.0
    for ptrange in list(sfs["NUM_MediumID_DEN_TrackerMuons"]["abseta_pt"][sel_eta].keys()):
        if pt >= float(ptrange.split("[")[1].split(",")[0]) and pt <= float(ptrange.split("]")[0].split(",")[1]):
            sel_pt = ptrange

    try:
        sf = sfs["NUM_MediumID_DEN_TrackerMuons"]["abseta_pt"][sel_eta][sel_pt]["value"] 
    except:
        #print("invalid low pt muon "," pt ",pt," eta ",eta)
        return 1.0
        

    if sf==0:
        sf=1.0
    #print("muon pt ",pt," low muon sf ",sf)
    return sf 

#for fake factor return the bin value from event content
#pt is a histogram
def ptFun(pt,numpyArr):
    newArr = np.full(len(numpyArr),1.0)
    bins = np.full(len(numpyArr),1)
    bins = np.vectorize(pt.FindBin)(numpyArr)
    bins = bins.astype(int)
    print(bins)
    for i,entry in enumerate(bins):
        newArr[i] = pt.GetBinContent(int(entry))
    #newArr = np.vectorize(pt.GetBinContent)(bins)

    return newArr

#def getPUfactor(masterArray):
#    import time
#    from ROOT import TFile, TTree, TH1D, TCanvas, TLorentzVector  
#    import numpy as np
#    tStart = time.time()
#    #print("Entering pileup function")
#    pileUpArray = np.full(len(masterArray['evt']),1.0)
#    #fData = ROOT.TFile.Open("${CMSSW_BASE}/src/ZH_Run2/pileup/data_pileup_{0:s}.root".format(args.year))
#    year = str(args.year)
#    if year=="2016":
#        fData = ROOT.TFile.Open("${CMSSW_BASE}/src/ZH_Run2/pileup/data_pileup_2016.root")
#    if year=="2017":
#        fData = ROOT.TFile.Open("${CMSSW_BASE}/src/ZH_Run2/pileup/data_pileup_2017.root")
#    if year=="2018":
#        fData = ROOT.TFile.Open("${CMSSW_BASE}/src/ZH_Run2/pileup/data_pileup_2018.root")
#
#    hData = fData.Get('pileup')
#    #print("hData={0:s}".format(str(hData)))
#    binWidth = hData.GetBinWidth(1)
#    xMin = hData.GetBinLowEdge(1)
#    nBins = hData.GetNbinsX()
#    xMax = xMin + nBins*binWidth
#    bins = np.linspace(xMin+0.5*binWidth,xMax-0.5*binWidth,nBins)
#    #print("nBins={0:d} binWidth={1:f} xMin={2:f} xMax={3:f}".format(nBins,binWidth,xMin,xMax))
#
#    hMC = TH1D("hMC","hMC",nBins,xMin,xMax)
#    hWeight = TH1D("hWeight","hWeight",nBins,xMin,xMax)
#    root_numpy.fill_hist(hMC,masterArray["nPU"])
#    root_numpy.fill_hist(hWeight,masterArray["nPU"],masterArray["Generator_weight"])
#    
#    nData = hData.GetSumOfWeights()
#    pData = np.array(hData)[1:-1]/nData
#    nMC = hMC.GetSumOfWeights()
#    pMC = np.array(hMC)[1:-1]
#    pMC /= nMC
#    pMC = np.maximum(1.e-5*np.ones_like(pMC),pMC)
#    #print("sum of pMC={0:f}".format(np.sum(pMC)))
#    ratio = np.divide(pData,pMC)
#    #print("ratio of data to mc numpy ",ratio[98])
#    minmask = np.where(masterArray["nPUtrue"]>98)
#    masterArray["nPUtrue"][minmask]=98
#    #print("master array after ceiling ",masterArray["nPUtrue"])
#    pileuptrue = masterArray["nPUtrue"]
#    for ipu,pu in enumerate(pileuptrue):
#        #print("pileup ",int(pileuptrue[ipu]))
#        #print("scale factor ",ratio[int(pileuptrue[ipu])])
#        pileUpArray[ipu] = ratio[int(pileuptrue[ipu])]
#    #print("After pileup loop:  time={0:.1f} s   time/event={1:.1f} us".format(time.time()-tStart,1.e6*(time.time()-tStart)))
#    return pileUpArray 


def getEventWeightDicitonary(year):
    #import copyreg, copy, pickle # for picking the bound methods due to the multiplrocess tool

    
    #from TauPOG.TauIDSFs.TauIDSFTool import TauIDSFTool
    #from TauPOG.TauIDSFs.TauIDSFTool import TauESTool
    #from TauPOG.TauIDSFs.TauIDSFTool import TauFESTool



    #maybe I need this per channel?

    #electrons
    eleES = cset_ele["UL-Electron-ID-SF"]
    #5 mmet
    #6 mmmt
    #7 mmtt
    #8 mmem

    EventWeights={
            "3_tauSF_7":[[["cat","==",7],["gen_match_3","==",5]],[np.vectorize(getAntiJetSF),["pt_3","decayMode_3"]]],
            "4_tauSF_7":[[["cat","==",7],["gen_match_4","==",5]],[np.vectorize(getAntiJetSF),["pt_4","decayMode_4"]]],
            "3_tauSF_6":[[["cat","==",6],["gen_match_3","==",5]],[np.vectorize(getAntiJetSF),["pt_3","decayMode_3"]]],
            "4_tauSF_6":[[["cat","==",6],["gen_match_4","==",5]],[np.vectorize(getAntiJetSF),["pt_4","decayMode_4"]]],
            "3_tauSF_5":[[["cat","==",5],["gen_match_3","==",5]],[np.vectorize(getAntiJetSF),["pt_3","decayMode_3"]]],
            "4_tauSF_5":[[["cat","==",5],["gen_match_4","==",5]],[np.vectorize(getAntiJetSF),["pt_4","decayMode_4"]]],


            "3_muonSF_8":[[["cat","==",8],["gen_match_3","==",15]],[np.vectorize(getAntiMuSF),["eta_3","gen_match_3"]]],
            "3_muonSF_7":[[["cat","==",7],["gen_match_3","==",15]],[np.vectorize(getAntiMuSF),["eta_3","gen_match_3"]]],
            "4_muonSF_7":[[["cat","==",7],["gen_match_4","==",15]],[np.vectorize(getAntiMuSF),["eta_4","gen_match_4"]]],
            "3_muonSF_6":[[["cat","==",6],["gen_match_3","==",15]],[np.vectorize(getAntiMuSF),["eta_3","gen_match_3"]]],
            "4_muonSF_6":[[["cat","==",6],["gen_match_4","==",15]],[np.vectorize(getAntiMuSF),["eta_4","gen_match_4"]]],
            "3_muonSF_5":[[["cat","==",5],["gen_match_3","==",15]],[np.vectorize(getAntiMuSF),["eta_3","gen_match_3"]]],
            "4_muonSF_5":[[["cat","==",5],["gen_match_4","==",15]],[np.vectorize(getAntiMuSF),["eta_4","gen_match_4"]]],

            "4_electronSF_8":[[["cat","==",8],["gen_match_4","==",15]],[np.vectorize(getAntiEleSF),["eta_4","gen_match_4"]]],
            "3_electronSF_7":[[["cat","==",7],["gen_match_3","==",15]],[np.vectorize(getAntiEleSF),["eta_3","gen_match_3"]]],
            "4_electronSF_7":[[["cat","==",7],["gen_match_4","==",15]],[np.vectorize(getAntiEleSF),["eta_4","gen_match_4"]]],
            "3_electronSF_6":[[["cat","==",6],["gen_match_3","==",15]],[np.vectorize(getAntiEleSF),["eta_3","gen_match_3"]]],
            "4_electronSF_6":[[["cat","==",6],["gen_match_4","==",15]],[np.vectorize(getAntiEleSF),["eta_4","gen_match_4"]]],
            "3_electronSF_5":[[["cat","==",5],["gen_match_3","==",15]],[np.vectorize(getAntiEleSF),["eta_3","gen_match_3"]]],
            "4_electronSF_5":[[["cat","==",5],["gen_match_4","==",15]],[np.vectorize(getAntiEleSF),["eta_4","gen_match_4"]]],


            #added later ... low muon pt scale factors from Muon Pog validated in Jsi region 
            #https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/Run2/UL
            "3_elecSF_8":[[["cat","==",8],["gen_match_3","==",15]],[np.vectorize(getElectronSF),["pt_3","eta_3"]]],
            "3_elecSF_5":[[["cat","==",5],["gen_match_3","==",15]],[np.vectorize(getElectronSF),["pt_3","eta_3"]]],
            "4_lowmuonSF_8":[[["cat","==",8],["gen_match_4","==",15]],[np.vectorize(getLowMuSF),["pt_4","eta_4"]]],
            "3_lowmuonSF_7":[[["cat","==",7],["gen_match_3","==",15]],[np.vectorize(getLowMuSF),["pt_3","eta_3"]]],
            "4_lowmuonSF_7":[[["cat","==",7],["gen_match_4","==",15]],[np.vectorize(getLowMuSF),["pt_4","eta_4"]]],
            "3_lowmuonSF_6":[[["cat","==",6],["gen_match_3","==",15]],[np.vectorize(getLowMuSF),["pt_3","eta_3"]]],
            "4_lowmuonSF_6":[[["cat","==",6],["gen_match_4","==",15]],[np.vectorize(getLowMuSF),["pt_4","eta_4"]]],
            "3_lowmuonSF_5":[[["cat","==",5],["gen_match_3","==",15]],[np.vectorize(getLowMuSF),["pt_3","eta_3"]]],
            "4_lowmuonSF_5":[[["cat","==",5],["gen_match_4","==",15]],[np.vectorize(getLowMuSF),["pt_4","eta_4"]]]
            }

    return EventWeights

def initialize(args):
    import os
    from copy import copy
    import yaml


    #Structure for mapping between processes and weights
    from utils.Weights import CommonWeights


    #sys.path.append('../SFs/')
    #import utils.SFs.ScaleFactor as SF

    import io
    import os
    from shutil import copyfile


    #gather functions for computing variables in the event loop
    from utils.functions import functs
    from ROOT import gInterpreter

    #Structure for plotting variables
    from utils.Parametrization import Category#, DylanCategory
    from utils.Parametrization import Process

    #Gather the analysis datasets and info
    allcats={}
    HAA_processes={}
    finalDistributions={}
    filelist = {}
    weightHistoDict={}
    jetWeightMultiplicity={}
    EventWeights={}
    datadrivenPackage={}
    sampleDict = {}
    #csvfile = "MCsamples_2016_v7_ZZ.csv"
    #csvfile = "MCsamples_2016_v7.csv"
    #categories = "cat_"+channel+"_2016.yaml"
    #processes = "processes_special_"+channel+".yaml"
    #categories = "cat_mmet_2016.yaml"
    #processes = "processes_special_mmet.yaml"
    csvfile = args.csvfile
    categories = args.categories
    processes = args.processes

    #dir = "/eos/home-s/shigginb/HAA_ntuples/2016_v7/"
    #dir = "/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov7_basic_10_6_4/src/2016_v7/"
    dir = args.dir



    for line in open(csvfile,'r').readlines() :
            # structure for csv to sampleDict conversion
            #[nickname]        = [category,xsec,numberOfEvents,finishedEvents,idk?,DASDataset]
            #print "reading line "
            #print line
            if "#" in line[0]: continue
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
    with io.open(categories,'r') as catsyam:
        categories = yaml.load(catsyam)

    #loading fake factor and data driven methods
    with io.open(processes,'r') as prosyam:
        processes_special = yaml.load(prosyam)

    allcats={}

    for category in categories:
        #print(category)
        #print(categories[category]['name'])
        #if categories[category]['name']==args.channel+"_inclusive":
        #if categories[category]['name']=="mmet_inclusive":
        #allows for anchors in the yaml file .. only load categories with names 
        try: categories[category]['name']
        except: continue
        tempcat = Category()
        tempcat.name=categories[category]['name']
        tempcat.cuts=categories[category]['cuts']
        tempcat.newvariables=categories[category]['newvariables']
        tempcat.vars=categories[category]['vars']
        #tempcat.systematics=categories[category]['systematics']
        allcats[tempcat.name]=tempcat
        #print allcats[tempcat.name].cuts
        # if args.extract and (args.channel=="mmmt" or args.channel=="mmet":
        #     allcats[tempcat.name].cuts["categorycuts"].append(["AMass","<=",120.0])
        # if args.extract and args.channel=="mmtt":
        #     allcats[tempcat.name].cuts["categorycuts"].append(["AMass","<=",130.0])
        # if args.extract and args.channel=="mmem":
        #     allcats[tempcat.name].cuts["categorycuts"].append(["AMass","<=",110.0])

    #allcats = {name: DylanCategory(cat) for name, cat in categories.items()}


    #for cat in allcats.keys():
        #print(allcats[cat].name)
        #print(allcats[cat].cuts)
        #print(allcats[cat].systematics)

    #loading standard processes
    HAA_processes={}

    #MC for background

    #HAA_process = {name: Process(dir, name, info) for name, info in sampleDict.items()}
    weightHistoDict = {}
    jetWeightMultiplicity = {}
    nevents = 250000

    if not args.onlydata:
        if not (args.datadrivenHAA or args.datadrivenSM):
            for sample in sampleDict.keys():
                temppro = Process()
                temppro.nickname=sample
                temppro.file=dir+sample+"_"+args.year+".root"

                frooin = ROOT.TFile.Open(temppro.file,"read")
                #temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue","genweight":"Generator_weight"}
                temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue","genweight":"weight"}
                temppro.cuts={sampleDict[sample][0]:""}
                if "ggTo2mu2tau" in sample:
                    #if args.year==2016: nevents = 250000
                    #if args.year==2017: nevents = 350000
                    #if args.year==2017: nevents = 500000
                    if args.extract:
                        #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001)} # SM Higgs xsec [pb] x BR Haa x 5 for DataMC control plots
                        #temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"Generator_weight"} # SM Higgs xsec [pb] x BR Haa x 5 for DataMC control plots
                        temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"weight"} # SM Higgs xsec [pb] x BR Haa x 5 for DataMC control plots
                        #temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"Generator_weight"}

                    else:
                        #temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(137.5*31.05*0.00005)} # worked before
                        #temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"Generator_weight"} # worked before
                        temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"weight"} # worked before
                        #temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"Generator_weight"}

                    #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001*5.0)} # SM Higgs xsec x BR Haa x 5 for DataMC control plots
                    #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.001* 5.0)} # SM Higgs xsec x BR Haa  for signal extraction data MC control plots(AN 17029)
                HAA_processes[temppro.nickname]=temppro

        if (args.datadrivenHAA or args.datadrivenSM):
            for sample in sampleDict.keys():
                temppro = Process()
                temppro.nickname=sample
                #temppro.file=sample+"_2016.root"
                #temppro.file=dir+sample+"_2016.root"
                temppro.file=dir+sample+"_"+args.year+".root"
                #temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue","genweight":"Generator_weight"}
                temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue","genweight":"weight"}
                if args.channel=="mmtt":
                    #truetau = [
                    #            [["OR"],
                    #            ["gen_match_3","==",5],
                    #            ["gen_match_4","==",5]
                    #            ]]
                    truetau = [["gen_match_3","==",5],
                                ["gen_match_4","==",5]
                                ]
                if args.channel=="mmem":
                    #truetau = [ [["OR"],
                    #            ["gen_match_3","==",15]
                    #            ],
                    #            [["OR"],
                    #            ["gen_match_4","==",15]
                    #            ] ]
                    truetau = [["gen_match_3","==",15],
                                ["gen_match_4","==",15]
                                ]
                if args.channel=="mmmt" or args.channel=="mmet":
                    #truetau = [[["OR"],
                    #            ["gen_match_3","==",15],
                    #            ["gen_match_4","==",5]
                    #            ]]
                    truetau = [["gen_match_3","==",15],
                               ["gen_match_4","==",5]
                                ]
                temppro.cuts={sampleDict[sample][0]:truetau} #ONLY SELECT PRMOPT FOR MC!
                if "ggTo2mu2tau" in sample:
                    #if args.year==2016: nevents = 250000
                    #if args.year==2017: nevents = 350000
                    #if args.year==2017: nevents = 500000
                    #for visualization!
                    if args.extract:
                        #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001)} # SM Higgs xsec x BR Haa x 5 for DataMC control plots
                        #temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"Generator_weight"} # SM Higgs xsec [pb] x BR Haa x 5 for DataMC control plots
                        temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"weight"} # SM Higgs xsec [pb] x BR Haa x 5 for DataMC control plots
                        #temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"Generator_weight"}

                    else:
                        #temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(137.5*31.05*0.00005)} # worked before
                        #temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"Generator_weight"} # worked before
                        temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"weight"} # worked before
                        #temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"theoryXsec":(48.37*0.001),"PU":"weightPUtrue","genweight":"Generator_weight"}
                HAA_processes[temppro.nickname]=temppro

        #OK make sure that process cuts are used in cutOnTreeArray procuts 
        if (args.datameasure):
            for proObj in HAA_processes.keys():
                if proObj!="data_obs" and proObj!="FF" and "ggTo2mu2tau" not in proObj:
                    if args.channel=="mmmt" or args.channel=="mmet":
                        #this create a new end distribution with like all the variables
                        if "ggTo2mu2tau" in proObj:
                            continue
                        HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","==",15]]
                        HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","==",5]]
                        HAA_processes[proObj].cuts["fake1"] = [["gen_match_3","!=",5]]
                        HAA_processes[proObj].cuts["fake2"] = [["gen_match_4","!=",5]]

                    if args.channel=="mmtt":
                        if "ggTo2mu2tau" in proObj:
                            continue
                        HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","==",5]]
                        HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","==",5]]
                        HAA_processes[proObj].cuts["fake1"] = [["gen_match_3","!=",5]]
                        HAA_processes[proObj].cuts["fake2"] = [["gen_match_4","!=",5]]
                    if args.channel=="mmem":
                        HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","!=",15]]
                        HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","!=",15]]
                        if "ggTo2mu2tau" in proObj:
                            continue
                        HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","==",15]]
                        HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","==",15]]
                    else: # default case is close enough
                        HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","==",0]]
                        HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","==",0]]
                        if "ggTo2mu2tau" in proObj:
                            continue
                        HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","!=",0]]
                        HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0]]
        elif (args.datadrivenHAA or args.datadrivenSM):
            for proObj in HAA_processes.keys():
                if proObj!="data_obs" and proObj!="FF":
                    #HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","==",0]]
                    #HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","==",0]]
                    #if "ggTo2mu2tau" in proObj:
                    #    print "skipping signal prompt mc"
                    #    continue  #this hasn't been added yet ... this is needed to omit signal!!
                    #HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","!=",0]]
                    #HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0]]
                    for key,subprocess in HAA_processes[proObj].cuts.items():
                        subprocess.append(["gen_match_3","!=",0])
                        subprocess.append(["gen_match_4","!=",0])


    #loading special processes ... fake factor and data!
    for process in processes_special:
        try: processes_special[process]['nickname']
        except: continue
        if(args.onlydata and processes_special[process]['nickname']!="data"):continue
        if(args.datameasure and processes_special[process]['nickname']!="data"):continue
        temppro = Process()
        temppro.nickname=processes_special[process]['nickname']
        temppro.cuts=processes_special[process]['cuts']
        temppro.weights=processes_special[process]['weights']
        temppro.file=dir+processes_special[process]['file']
        HAA_processes[temppro.nickname]=temppro
    #print "==========================================================================================="
    #print " All the processes !!!!!  ",HAA_processes
    #print "==========================================================================================="


    if not args.onlydata:
        #precise weights for jet multiplicity
        from utils.Weights import jet_exc_samples
        from utils.Weights import jetIncOnly
        inc_samps = jetIncOnly[args.year]
        jet_exc_samples = jet_exc_samples[args.year]

        #file list for easy reading
        for proObj in HAA_processes.keys():
            #if "2016" in args.year:
            #    filelist[proObj+"_pre"]=HAA_processes[proObj].file
            #    filelist[proObj+"_post"]=HAA_processes[proObj].file.replace("pre","post")
            #else:
            filelist[proObj]=HAA_processes[proObj].file
        #if "2016" in args.year: print("ok... adding pre and post for combined SoW  \n",filelist)
                

        # adding kfactor to relevent samples
        for sample in sampleDict.keys():
            if sample in jet_exc_samples and sample.startswith("DY"):
                HAA_processes[sample].weights.update({"kfactor":1.1637})
            if sample in jet_exc_samples and sample.startswith("W"):
                HAA_processes[sample].weights.update({"kfactor":1.221})

        weightHistoDict = {}
        for nickname, filename in filelist.items():
            #frooin = ROOT.TFile(dir+filename)
            frooin = ROOT.TFile(filename)
            hists = None
            try:
                #if nickname.startswith("DYJetsToLL"):
                #    hists=[frooin.Get("hWeights").Clone(),
                #           frooin.Get("DY1genWeights").Clone(),
                #           frooin.Get("DY2genWeights").Clone(),
                #           frooin.Get("DY3genWeights").Clone(),
                #           frooin.Get("DY4genWeights").Clone(),
                #           ]
                #elif nickname.startswith("WJetsToLNu"):
                #    hists=[frooin.Get("hWeights").Clone(),
                #          frooin.Get("W1genWeights").Clone(),
                #          frooin.Get("W2genWeights").Clone(),
                #          frooin.Get("W3genWeights").Clone(),
                #          frooin.Get("W4genWeights").Clone(),
                #          ]
                #else:
                hists = frooin.Get("hWeights").Clone()
                hists.SetDirectory(0)
                if isinstance(hists, list):
                    for hist in hists:
                        hist.SetDirectory(0)
                weightHistoDict[nickname] = hists
                frooin.Close()
            except:
                print(nickname," hWeights prob doesn't exist?")

            if "2016" in args.year: 
                if "pre" in args.year:
                    frooin = ROOT.TFile(filename.replace("pre","post"))
                if "post" in args.year:
                    frooin = ROOT.TFile(filename.replace("post","pre"))
                hists = None
                try:
                    #if nickname.startswith("DYJetsToLL"):
                    #    hists=[frooin.Get("hWeights").Clone(),
                    #           frooin.Get("DY1genWeights").Clone(),
                    #           frooin.Get("DY2genWeights").Clone(),
                    #           frooin.Get("DY3genWeights").Clone(),
                    #           frooin.Get("DY4genWeights").Clone(),
                    #           ]
                    #elif nickname.startswith("WJetsToLNu"):
                    #    hists=[frooin.Get("hWeights").Clone(),
                    #          frooin.Get("W1genWeights").Clone(),
                    #          frooin.Get("W2genWeights").Clone(),
                    #          frooin.Get("W3genWeights").Clone(),
                    #          frooin.Get("W4genWeights").Clone(),
                    #          ]
                    #else:
                    hists = frooin.Get("hWeights").Clone()
                    hists.SetDirectory(0)
                    if isinstance(hists, list):
                        for hist in hists:
                            hist.SetDirectory(0)
                    if "pre" in args.year:
                        weightHistoDict[nickname+"post"] = hists
                    if "post" in args.year:
                        weightHistoDict[nickname+"pre"] = hists
                    frooin.Close()
                except:
                    print(nickname+"_post"," hWeights prob doesn't exist?")

        #print weightHistoDict
        #for key in weightHistoDict.keys():
            #if weightHistoDict[key].Get
        #print weightHistoDict["DYJetsToLLext1"]
        #print weightHistoDict["WJetsToLNu"]


        jetWeightMultiplicity = {}
        NJetWeights = {}
        # for all the DYJet and WJet files make sure that the exts are summed PRIOR to makeCutsOnTreeArray
        #print inc_samps
        #in ultra legacy only focusing on datadriven so I am not using ALL the DY and WJets files anymore
        #cancelDYandWxsecreweight = 0
        cancelDYandWxsecreweight = 1
        #for inc_group in inc_samps:
        #    sumOfWeights = np.zeros(5)
        #    for sample in inc_group:
        #        try:
        #            weights = weightHistoDict[sample]
        #        except:
        #            cancelDYandWxsecreweight = 1
        #            break
        #        sumOfWeights += [w.GetSumOfWeights() for w in weights]
        #    if cancelDYandWxsecreweight: break
        #    sumOfWeights[sumOfWeights==0] = 1e-10 # if value is 0 to stop infs
        #    #print(sumOfWeights)
        #    if not cancelDYandWxsecreweight:
        #        for sample in inc_group:
        #            jetWeightMultiplicity[sample] = sumOfWeights


        #if args.year=="2018" and not cancelDYandWxsecreweight:
        #    DY1JetsFile = ROOT.TFile.Open(filelist["DY1JetsToLL"],"read")
        #    jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
        #    DY2JetsFile = ROOT.TFile.Open(filelist["DY2JetsToLL"],"read")
        #    jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
        #    DY3JetsFile = ROOT.TFile.Open(filelist["DY3JetsToLL"],"read")
        #    jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
        #    DY4JetsFile = ROOT.TFile.Open(filelist["DY4JetsToLL"],"read")
        #    jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()

        #    W1JetsFile = ROOT.TFile.Open(filelist["W1JetsToLNu"],"read")
        #    jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()

        #    W2JetsFile = ROOT.TFile.Open(filelist["W2JetsToLNu"],"read")
        #    jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()

        #    W3JetsFile = ROOT.TFile.Open(filelist["W3JetsToLNu"],"read")
        #    jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()

        #    W4JetsFile = ROOT.TFile.Open(filelist["W4JetsToLNu"],"read")
        #    jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()

        #if args.year=="2017" and not cancelDYandWxsecreweight:
        #    DY1JetsFile = ROOT.TFile.Open(filelist["DY1JetsToLL"],"read")
        #    jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
        #    DY1JetsFile_ext1 = ROOT.TFile.Open(filelist["DY1JetsToLL_ext1"],"read")
        #    jetWeightMultiplicity["DY1JetsToLL_ext1"]=DY1JetsFile_ext1.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["DY1JetsToLL_ext1"]+=jetWeightMultiplicity["DY1JetsToLL"]
        #    jetWeightMultiplicity["DY1JetsToLL"]+=jetWeightMultiplicity["DY1JetsToLL_ext1"]

        #    DY2JetsFile = ROOT.TFile.Open(filelist["DY2JetsToLL"],"read")
        #    jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
        #    DY2JetsFile_ext1 = ROOT.TFile.Open(filelist["DY2JetsToLL_ext1"],"read")
        #    jetWeightMultiplicity["DY2JetsToLL_ext1"]=DY2JetsFile_ext1.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["DY2JetsToLL_ext1"]+=jetWeightMultiplicity["DY2JetsToLL"]
        #    jetWeightMultiplicity["DY2JetsToLL"]+=jetWeightMultiplicity["DY2JetsToLL_ext1"]

        #    DY3JetsFile = ROOT.TFile.Open(filelist["DY3JetsToLL"],"read")
        #    jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
        #    DY3JetsFile_ext1 = ROOT.TFile.Open(filelist["DY3JetsToLL_ext1"],"read")
        #    jetWeightMultiplicity["DY3JetsToLL_ext1"]=DY3JetsFile_ext1.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["DY3JetsToLL_ext1"]+=jetWeightMultiplicity["DY3JetsToLL"]
        #    jetWeightMultiplicity["DY3JetsToLL"]+=jetWeightMultiplicity["DY3JetsToLL_ext1"]

        #    DY4JetsFile = ROOT.TFile.Open(filelist["DY4JetsToLL"],"read")
        #    jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()

        #    W1JetsFile = ROOT.TFile.Open(filelist["W1JetsToLNu"],"read")
        #    jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()

        #    W2JetsFile = ROOT.TFile.Open(filelist["W2JetsToLNu"],"read")
        #    jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()

        #    W3JetsFile = ROOT.TFile.Open(filelist["W3JetsToLNu"],"read")
        #    jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()

        #    W4JetsFile = ROOT.TFile.Open(filelist["W4JetsToLNu"],"read")
        #    jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()

        #if "2016" in args.year and not cancelDYandWxsecreweight:
        #    DY1JetsFile = ROOT.TFile.Open(filelist["DY1JetsToLL"],"read")
        #    jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
        #    DY2JetsFile = ROOT.TFile.Open(filelist["DY2JetsToLL"],"read")
        #    jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
        #    DY3JetsFile = ROOT.TFile.Open(filelist["DY3JetsToLL"],"read")
        #    jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
        #    DY4JetsFile = ROOT.TFile.Open(filelist["DY4JetsToLL"],"read")
        #    jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()

        #    W1JetsFile = ROOT.TFile.Open(filelist["W1JetsToLNu"],"read")
        #    jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()

        #    W2JetsFile = ROOT.TFile.Open(filelist["W2JetsToLNu"],"read")
        #    jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()
        #    W2JetsFileext1 = ROOT.TFile.Open(filelist["W2JetsToLNu_ext1"],"read")
        #    jetWeightMultiplicity["W2JetsToLNu"]+=W2JetsFileext1.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["W2JetsToLNu_ext1"]=W2JetsFileext1.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["W2JetsToLNu_ext1"]+=W2JetsFile.Get("hWeights").GetSumOfWeights()

        #    W3JetsFile = ROOT.TFile.Open(filelist["W3JetsToLNu"],"read")
        #    jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()
        #    W3JetsFileext1 = ROOT.TFile.Open(filelist["W3JetsToLNu_ext1"],"read")
        #    jetWeightMultiplicity["W3JetsToLNu"]+=W3JetsFileext1.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["W3JetsToLNu_ext1"]=W3JetsFileext1.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["W3JetsToLNu_ext1"]+=W3JetsFile.Get("hWeights").GetSumOfWeights()

        #    W4JetsFile = ROOT.TFile.Open(filelist["W4JetsToLNu"],"read")
        #    W4JetsFileext1 = ROOT.TFile.Open(filelist["W4JetsToLNu_ext1"],"read")
        #    W4JetsFileext2 = ROOT.TFile.Open(filelist["W4JetsToLNu_ext2"],"read")
        #    jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["W4JetsToLNu"]+=W4JetsFileext1.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["W4JetsToLNu"]+=W4JetsFileext2.Get("hWeights").GetSumOfWeights()

        #    jetWeightMultiplicity["W4JetsToLNu_ext1"]=W4JetsFileext1.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["W4JetsToLNu_ext1"]+=W4JetsFile.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["W4JetsToLNu_ext1"]+=W4JetsFileext2.Get("hWeights").GetSumOfWeights()

        #    jetWeightMultiplicity["W4JetsToLNu_ext2"]=W4JetsFileext2.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["W4JetsToLNu_ext2"]+=W4JetsFile.Get("hWeights").GetSumOfWeights()
        #    jetWeightMultiplicity["W4JetsToLNu_ext2"]+=W4JetsFileext1.Get("hWeights").GetSumOfWeights()


        Bkg = ["DY","W","TT","ST","EWK"]
        #irBkg = ["ZZ","ZHToTauTau","vbf","WHTT"]
        irBkg = ["ZZ"]
        #TrialphaBkg = ["ttZ","ttW","WWZ","WZZ","ZZZ","WWW_4F","HZJ"]
        TrialphaBkg = ["ttZ","WWZ","WZZ","ZZZ","HZJ"]
        #rareBkg = ["Other","rare","WZ"]
        rareBkg = ["Other"]
        finalDistributions = {}
        finalDistributions["Bkg"]=Bkg
        finalDistributions["data_obs"]=["data_obs"]
        finalDistributions["SS_relaxed_data"]=["SS_relaxed_data"]
        finalDistributions["a15"]=["a15"]
        finalDistributions["a20"]=["a20"]
        finalDistributions["a25"]=["a25"]
        finalDistributions["a30"]=["a30"]
        finalDistributions["a35"]=["a35"]
        finalDistributions["a40"]=["a40"]
        finalDistributions["a45"]=["a45"]
        finalDistributions["a50"]=["a50"]
        finalDistributions["a55"]=["a55"]
        finalDistributions["a60"]=["a60"]
        finalDistributions["irBkg"]=irBkg
        finalDistributions["TrialphaBkg"]=TrialphaBkg
        finalDistributions["rareBkg"]=rareBkg
        if args.datameasure:
            finalDistributions["prompt1"]=["prompt1"]
            finalDistributions["prompt2"]=["prompt2"]
        if args.datadrivenHAA or args.datadrivenSM:
            Bkg = ["FF"]
            #rareBkg = ["DY","W","TT","ST","EWK","Other","rare","WZ"]
            finalDistributions["Bkg"]=Bkg
            #finalDistributions["rareBkg"]=rareBkg
    else: #only data
        finalDistributions = {}
        finalDistributions["data_obs"]=["data_obs"]
        Bkg = ["FF"]
        finalDistributions["Bkg"]=Bkg

    allSkims = {}
    finalSkims ={}

    datadrivenPackage={}
    datadrivenPackage["bool"]=False

    fakefactorObj= fakefactor()
    if args.datadrivenSM:
        #fakefactorObj = fakefactor("/eos/home-s/shigginb/fakefactor/")
        #fakefactorObj.loadHistograms("/eos/home-s/shigginb/fakefactors/",args.year)
        fakefactorObj.loadHistograms("./fakefactors/",args.year)
        #fakefactorObj.Print()
    if args.datadrivenHAA:
        datadrivenPackage={}
        datadrivenPackage["bool"]=args.datadrivenHAA
        fakemeasurefile_1 = ROOT.TFile.Open(str(args.ffin1),"READ")
        fakemeasurefile_1.cd()
        datadrivenPackage["fakerate1"]=fakemeasurefile_1.Get("fakerate1")
        datadrivenPackage["fitrate1"]=fakemeasurefile_1.Get("tf_1")

        fakemeasurefile_2 = ROOT.TFile.Open(str(args.ffin2),"READ")
        fakemeasurefile_2.cd()
        datadrivenPackage["fakerate2"]=fakemeasurefile_2.Get("fakerate1") 
        datadrivenPackage["fitrate2"]=fakemeasurefile_2.Get("tf_1")
        datadrivenPackage["file1"] = fakemeasurefile_1
        datadrivenPackage["file2"] = fakemeasurefile_2

    #exit()
    # ROOT.fail
    return allcats, HAA_processes,finalDistributions,weightHistoDict,jetWeightMultiplicity,datadrivenPackage,fakefactorObj





def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

def fstar(args):
    #print("arggggg matey!")
    return f(*args)

def slimskimstar(args):
    #return slimskim(*args)
    print("in slimskimstar")
    return slimskimoutput(*args)

def f(process, categories):
    info('function f')
    #print('hello', process.nickname)
    #for cat in categories.keys():
        #print("category ",cat)
    return 1

def createFakeFactorHistos(allcats, inputFFile,promptsub):
    newVarVals={}
    if promptsub:
        treetypes = ["Nominal_data_obs","Nominal_prompt1","Nominal_prompt2"]
    else:
        treetypes = ["Nominal_data_obs"]
    # for cat in allcats.keys():
    #     if not "_inclusive" in cat:
    #         histodict={}
    #         for treename in treetypes:
    #             histodict[cat]={}
    #             histodict[cat][treename]={}
    histodict = {c:{t: dict() for t in treetypes} for c in allcats.keys()}
    #print "histodict ",histodict


    #creating the Fake Factor histograms from pre-defined numpy arrays
    for cat in histodict.keys():
        #print cat
        if "inclusive" in cat: continue
        for treename in treetypes:
            tree = inputFFile[cat][treename]
            fakefactorArray = tree.arrays(library="np")
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





def makeCutsOnTreeArray(processObj, inputArray,allcats,weightHistoDict,systematic,jetWeightMultiplicity,fakefactorObj,args):

    from utils.functions import functs
    from utils.Weights import CommonWeights
    from utils.Weights import jet_exc_samples
    import pickle

    jet_exc_samples = jet_exc_samples[args.year]

    from ROOT import gInterpreter
    import copy
    commonweight = CommonWeights["lumi"+args.year][0]
    print("working on process obj ",processObj.nickname)

    for process in processObj.cuts.keys():
        procuts = processObj.cuts[process]
        #print "working on process ",process
        for cat in allcats.keys():
            skimArray = {}
            #masterArray = inputArray.copy()
            masterArray = inputArray.copy()
            cuts=[]
            #print "working on category ",cat
            #print "with vars ",allcats[cat].vars
            #print "starting length of dictionary ",len(masterArray["mll"])

            for cuttype in allcats[cat].cuts.keys():
                for cut in allcats[cat].cuts[cuttype]:
                    cuts.append(cut)
            if procuts!="":
                for pc in procuts:
                    cuts.append(pc)
            if args.extract and args.channel=="mmmt":
                cuts.append(["AMass","<=",120.0])
                cuts.append(["mll-mtt",">",0.0])
            elif args.extract and args.channel=="mmtt":
                cuts.append(["AMass","<=",135.0])
                cuts.append(["mll-mtt",">",0.0])
                #if "2016" in args.year:
                #else:
                #    cuts.append(["AMass","<=",130.0])
                #    cuts.append(["mll-mtt",">",0.0])
            elif args.extract and args.channel=="mmet":
                cuts.append(["AMass","<=",120.0])
                cuts.append(["mll-mtt",">",0.0])
            elif args.extract and args.channel=="mmem":
                cuts.append(["AMass","<=",110.0])
                cuts.append(["mll-mtt",">",0.0])
            #print(process)
            #print(process.cuts)
            #print(masterArray.keys())
            plottedVars = []
            newVarVals={}
            for variableHandle in allcats[cat].vars.keys():
                variable = allcats[cat].vars[variableHandle][0]
                if "[" in variable:
                    #print "adding jagged variable to array ",variable
                    basevar = variable.split("[")[0]
                    index = int(variable.split("[")[1].split("]")[0])
                    val = masterArray[basevar][:,index]
                    #masterArray[basevar+"_"+str(index)]=val
                    #plottedVars.append(basevar+"_"+str(index))
                    masterArray[variableHandle+"_"+str(index)]=val
                    plottedVars.append(variableHandle+"_"+str(index))
                    plottedVars.append(variable)
                else:
                    plottedVars.append(variableHandle)
                    plottedVars.append(variable)

            for var in allcats[cat].newvariables.keys():
                newVarVals[var]=0.0
                plottedVars.append(var) # new vars just label
            for var in newVarVals.keys():
                arguments = allcats[cat].newvariables[var][2]
                tempvals=[]
                for ag in arguments:
                    tempvals.append(returnArray(masterArray,ag))
                masterArray[var] = functs[allcats[cat].newvariables[var][0]](*tempvals)


            plottedVars.append("nickname")
            #masterArray["nickname"]=np.full(len(masterArray['evt']),str(processObj.nickname))
            masterArray["nickname"]=np.full(len(masterArray['mll']),str(processObj.nickname))
            masterArray["nickname"]=masterArray["nickname"].astype("S40")
            #print masterArray["nickname"]

            if (process=="data_obs"):
                masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)

                print(cuts)
                mask = cutOnArray(masterArray,cuts)
                #masterArray["mask"]=mask
                masterArray["finalweight"] *= mask.astype(int)
                weightfinal = 1.0   #don't weight the data!!

                skipEvents = np.where(mask==0)[0]
                skimArray={}
                #print("before skim", len(masterArray["finalweight"]))
                for key in masterArray.keys():
                    try:
                        skimArray[key] = masterArray[key][mask]
                    except:
                        #print "length problem? length of key in master ",len(masterArray[key])," length of mask ",len(mask)
                        #print "skipping branch ",key
                        continue
                #print("after skim", len(skimArray["mll"]), processObj.file)
                print("cat ",cat,"  events ",len(skimArray["mll"]))
                #if len(skimArray["mll"])==0:
                #    continue
                print("data obs","  ",cat," events ",len(masterArray["finalweight"])," data finalyield ",np.sum(masterArray["finalweight"]))

                #for key in skimArray.keys():
                for key in list(skimArray):
                    if key not in plottedVars and key != "finalweight" and key != "evt":
                        del skimArray[key]
                #skimArrayPerCat[systematic+":"+cat+":"+processObj.nickname+":"+process] = skimArray

            if (process=="SS_relaxed_data"):
                masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)
                #print("events before cuts",len(masterArray["evt"]))
                cuts = HAA_processes["SS_relaxed_data"].cuts["SS_relaxed_data"]
                mask = cutOnArray(masterArray,cuts)
                #masterArray["mask"]=mask
                masterArray["finalweight"] *= mask.astype(int)
                weightfinal = 1.0   #don't weight the data!!

                skipEvents = np.where(mask==0)[0]
                skimArray={}
                for key in masterArray.keys():
                    try:
                        skimArray[key] = masterArray[key][mask]
                    except:
                        continue
                if len(skimArray["mll"])==0:
                    continue

                #for key in skimArray.keys():
                for key in list(skimArray):
                    if key not in plottedVars and key != "finalweight" and key != "evt":
                        del skimArray[key]
                #print("events after cuts",len(skimArray["mll"]))
                print("SS relaxed ","  ",cat," events ",len(masterArray["finalweight"])," data finalyield ",np.sum(masterArray["finalweight"]))
                #skimArrayPerCat[systematic+":"+cat+":"+processObj.nickname+":"+process] = skimArray

            

            if(process=="FF" and datadrivenPackage["bool"] and args.datadrivenHAA):
                print("running FF")
                masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)


                tempmask=np.full(len(masterArray["evt"]),1.0)

                #the actual events that pass the FF_1 criteria
                cuts_1 = HAA_processes["FF"].cuts["FF_1"]
                cuts_2 = HAA_processes["FF"].cuts["FF_2"]
                cuts_12 = HAA_processes["FF"].cuts["FF_12"]
                if args.extract and (args.channel=="mmmt" or args.channel=="mmet"):
                    cuts_1.append(["AMass","<=",120.0])
                    cuts_1.append(["mll-mtt",">",0.0])
                    cuts_2.append(["AMass","<=",120.0])
                    cuts_2.append(["mll-mtt",">",0.0])
                    cuts_12.append(["AMass","<=",120.0])
                    cuts_12.append(["mll-mtt",">",0.0])
                if args.extract and args.channel=="mmtt":
                    if "2016" in args.year:
                        cuts_1.append(["AMass","<=",135.0])
                        cuts_1.append(["mll-mtt",">",0.0])
                        cuts_2.append(["AMass","<=",135.0])
                        cuts_2.append(["mll-mtt",">",0.0])
                        cuts_12.append(["AMass","<=",135.0])
                        cuts_12.append(["mll-mtt",">",0.0])
                    else:
                        cuts_1.append(["AMass","<=",130.0])
                        cuts_1.append(["mll-mtt",">",0.0])
                        cuts_2.append(["AMass","<=",130.0])
                        cuts_2.append(["mll-mtt",">",0.0])
                        cuts_12.append(["AMass","<=",130.0])
                        cuts_12.append(["mll-mtt",">",0.0])
                if args.extract and args.channel=="mmem":
                    cuts_1.append(["AMass","<=",110.0])
                    cuts_1.append(["mll-mtt",">",0.0])
                    cuts_2.append(["AMass","<=",110.0])
                    cuts_2.append(["mll-mtt",">",0.0])
                    cuts_12.append(["AMass","<=",110.0])
                    cuts_12.append(["mll-mtt",">",0.0])
                tempmask_1 = cutOnArray(masterArray,cuts_1)
                #print tempmask_1[:1000]
                tempmask_2 = cutOnArray(masterArray,cuts_2)
                #print tempmask_2[:1000]
                tempmask_12 = cutOnArray(masterArray,cuts_12)
                #print tempmask_12[:1000]


                #FF_1
                #causal catching ... the pt may be outside the shape from the histogram... if so we need the constant fit value for extrapolation
                print("using fake rates ",datadrivenPackage.keys())
                fitmask_1 = cutOnArray(masterArray,[["pt_3","<",datadrivenPackage["fakerate1"].GetBinLowEdge(datadrivenPackage["fakerate1"].GetNbinsX())],["pt_3",">",datadrivenPackage["fakerate1"].GetBinLowEdge(2)]])
                fitmask_1 = fitmask_1.astype(int)

                ptarr_1 = masterArray["pt_3"]
                ffweight_1 = ptFun(datadrivenPackage["fakerate1"],ptarr_1)
                ffweight_1[fitmask_1]=datadrivenPackage["fitrate1"].GetParameter(0)
                ffweight_1 = ffweight_1/(1.0000000001 - ffweight_1)

                #by single value
                #ffweight_1 = np.ones(len(tempmask_1))*(0.03/(1-0.03))
                #ff measured as T/L(including tight)
                #print("loading fit parms")
                #ffweight_1 = np.ones(len(tempmask_1))*(datadrivenPackage["fitrate1"].GetParameter(0)/(1-datadrivenPackage["fitrate1"].GetParameter(0)))

                #ff measured as T/L(excluding tight)
                #ffweight_1 = np.ones(len(tempmask_1))*datadrivenPackage["fitrate1"].GetParameter(0)


                #FF_2
                fitmask_2 = cutOnArray(masterArray,[["pt_4","<",datadrivenPackage["fakerate2"].GetBinLowEdge(datadrivenPackage["fakerate2"].GetNbinsX())],["pt_4",">",datadrivenPackage["fakerate2"].GetBinLowEdge(2)]])
                fitmask_2 = fitmask_2.astype(int)

                ptarr_2 = masterArray["pt_4"]
                ffweight_2 = ptFun(datadrivenPackage["fakerate2"],ptarr_2)
                ffweight_2[fitmask_2]=datadrivenPackage["fitrate2"].GetParameter(0)
                ffweight_2 = ffweight_2/(1.0000000001 - ffweight_2)

                ffweight = -1.0 * ffweight_1 * ffweight_2
                ffweight *= tempmask_12

                ffweight_1 *= tempmask_1
                ffweight_2 *= tempmask_2
                finalWeight = ffweight_1 + ffweight_2 + ffweight


                masterArray["finalweight"] *= finalWeight
                print("f1 contribution ",np.sum(tempmask_1),np.sum(ffweight_1))
                print("f2 contribution ",np.sum(tempmask_2),np.sum(ffweight_2))
                print("f12 contribution ",np.sum(tempmask_12),np.sum(ffweight))
                print("summed final weight ",np.sum(finalWeight))

                keepEvents = np.where(finalWeight!=0.0)[0]

                skimArray={}
                for key in masterArray.keys():
                    skimArray[key] = masterArray[key][keepEvents]

                #print("after skim", len(skimArray["mll"]), processObj.file)
                if len(skimArray["mll"])==0:
                    continue

                #for key in skimArray.keys():
                for key in list(skimArray):
                    if key not in plottedVars and key != "finalweight" and key != "evt":
                        del skimArray[key]
                #return skimArray
                print("FF","  ",cat," events ",len(masterArray["finalweight"])," data finalyield ",np.sum(masterArray["finalweight"]))
                #skimArrayPerCat[systematic+":"+cat+":"+processObj.nickname+":"+process] = skimArray

            if process not in ["data_obs","FF","FF_1","FF_2","FF_12","SS_relaxed_data"]:
                EventWeights = getEventWeightDicitonary(args.year)

                #masterArray['finalweight']=np.full(len(masterArray['evt']),1.0)
                masterArray['finalweight']=np.full(len(masterArray['mll']),1.0)


                mask = cutOnArray(masterArray,cuts)
                #masterArray["mask"]=mask
                masterArray["finalweight"] *= mask.astype(int)

                print("MC before cuts ",len(masterArray["finalweight"]))
                #print("length of mask ",len(mask))
                #print("length of array ",len(masterArray["finalweight"]))
                #for key,value in masterArray.items():
                    #masterArray[key] = value[mask]
                skimArray={}
                for key in masterArray.keys():
                    try:
                        #skimArray[key] = masterArray[key][mask]
                        masterArray[key] = masterArray[key][mask]
                    except:
                        print(key)
                print("MC after cuts ",len(masterArray["finalweight"]))

                weightfinal = 1.0   #don't weight the data!!

                skipEvents = np.where(mask==0)[0]
                if len(masterArray["mll"])==0:
                    continue

                weightDict = processObj.weights
                weightfinal = commonweight
                nickname = processObj.nickname
                for scalefactor in weightDict.keys():
                    if scalefactor == "kfactor":
                        weightfinal =  weightfinal * (1 / float(weightDict[scalefactor]))
                    elif scalefactor in ["PU","genweight"]:
                    #elif scalefactor in ["genweight"]:
                        masterArray["finalweight"] *= (returnArray(masterArray,weightDict[scalefactor]))
                    #elif scalefactor == "PU": # calculating the pileup myself
                    #    masterArray["finalweight"] *= getPUfactor(masterArray)
                    elif scalefactor =="theoryXsec":
                        weightfinal =  weightfinal * float(weightDict[scalefactor])
                print(scalefactor,"   ",sum(masterArray["finalweight"]))

                print("finalweight after PU and kFactor ",masterArray["finalweight"][:100])
                sumOfWeights=0.0
                #if nickname in jet_exc_samples:
                #    jetweights = 5*[0]
                #    if nickname.startswith("DY"):
                #        jet0Xsec = 4673.65
                #        if "2016" in args.year:
                #            jetweights[0] = jet0Xsec/jetWeightMultiplicity["DYJetsToLLext1"][0]
                #            for nj in range(1, 5):
                #                njname = "DY%dJetsToLL"%nj
                #                norm1 = jetWeightMultiplicity["DYJetsToLLext1"][nj]
                #                norm2 = jetWeightMultiplicity[njname]
                #                jetweights[nj] = float(weightDict["kfactor"])*HAA_processes[njname].weights["xsec"] / (norm1 + norm2)
                #        elif args.year=="2017":
                #            jetweights[0] = jet0Xsec/jetWeightMultiplicity["DYJetsToLL"][0]
                #            for nj in range(1, 5):
                #                njname = "DY%dJetsToLL"%nj
                #                #norm1 = jetWeightMultiplicity["DYJetsToLL_ext1"][nj]
                #                norm1 = jetWeightMultiplicity["DYJetsToLL"][nj]
                #                norm2 = jetWeightMultiplicity[njname]
                #                jetweights[nj] = float(weightDict["kfactor"])*HAA_processes[njname].weights["xsec"] / (norm1 + norm2)
                #        elif args.year=="2018":
                #            jetweights[0] = jet0Xsec/jetWeightMultiplicity["DYJetsToLL"][0]
                #            for nj in range(1, 5):
                #                njname = "DY%dJetsToLL"%nj
                #                norm1 = jetWeightMultiplicity["DYJetsToLL"][nj]
                #                norm2 = jetWeightMultiplicity[njname]
                #                jetweights[nj] = float(weightDict["kfactor"])*HAA_processes[njname].weights["xsec"] / (norm1 + norm2)
                #    else:
                #        jet0Xsec = 49033.2
                #        jetweights[0] = jet0Xsec/jetWeightMultiplicity["WJetsToLNu"][0]
                #        for nj in range(1, 5):
                #            njname = "W%dJetsToLNu"%nj
                #            norm1 = jetWeightMultiplicity["WJetsToLNu"][nj]
                #            norm2 = jetWeightMultiplicity[njname]
                #            jetweights[nj] = float(weightDict["kfactor"])*HAA_processes[njname].weights["xsec"] / (norm1 + norm2)

                #    for i_jet,weight in enumerate(jetweights):
                #        if args.year==2016:
                #            njetmask = masterArray["LHE_Njets"]==i_jet
                #        else:
                #            njetmask = masterArray["njets"]==i_jet
                #        masterArray["finalweight"] [njetmask] *= weight
                #        #print "events that pass ",i_jet," jets ",np.count_nonzero(masterArray["finalweight"] [njetmask])

                #    #print "jet weight array ",jetweights
                #    #print  " sample name ",nickname,"xsec ",HAA_processes[nickname].weights["xsec"]," events that pass ", np.count_nonzero(masterArray["finalweight"])
                if not type(weightHistoDict[nickname])==list:
                #if not type(weightHistoDict[nickname])==list:
                    sumOfWeights = 0.0

                    for nic,sowhist in weightHistoDict.items():
                        extensions = ["_ext1", "_ext2", "_ext3", "ext1", "ext2","pre","post"]
                        for ext in extensions:
                         #sum sow from all versions of this process.
                            if nickname.replace(ext,"")==nic or nic.replace(ext,"")==nickname:
                                #print "adding ",nickname,"   to   ",nic
                                sumOfWeights += sowhist.GetSumOfWeights()
                                break

                    # for nic,sowhist in weightHistoDict.items():
                    #     if nickname in nic:
                    #         sumOfWeights += sowhist.GetSumOfWeights()

                    if sumOfWeights != 0.0:
                        weightfinal = weightfinal * HAA_processes[nickname].weights["xsec"]/ sumOfWeights

                masterArray["finalweight"] *= weightfinal

                print(" sample name ",nickname,"xsec ",HAA_processes[nickname].weights["xsec"]," SoW ",sumOfWeights," events that pass ", np.count_nonzero(masterArray["finalweight"]))

                #multiply by scalar weight
                #print "finalweight before per event scaling ",masterArray["finalweight"][:100]


                #eventWeightDict = process.eventWeights
                eventWeightDict = EventWeights
                if eventWeightDict:
                    for scalefactor in eventWeightDict.keys():
                        #print("working on event sf ",scalefactor)

                        cutlist = eventWeightDict[scalefactor][0]
                        weightMask=cutOnArray(masterArray,cutlist)
                        weightMask=weightMask.astype(float)

                        if type(eventWeightDict[scalefactor][1][0])==float:
                            weightMask *= eventWeightDict[scalefactor][1][0]

                        if hasattr(eventWeightDict[scalefactor][1][0],'__call__'):
                            arguments = eventWeightDict[scalefactor][1][1]
                            tempvals=[]

                            for ag in arguments:
                                #print("input for event weight ",ag)
                                tempvals.append(returnArray(masterArray,ag))
                            weightMask*= eventWeightDict[scalefactor][1][0](*tempvals)

                        if scalefactor!="fake" and scalefactor!="fake1" and scalefactor!="fake2":
                            weightMask[np.where(weightMask==0.0)]=1.0


                        masterArray["finalweight"] *= weightMask
                        #print(scalefactor,"   ",sum(masterArray["finalweight"]))

                #skimArray["finalweight"] = masterArray["finalweight"][mask]
                #skimArray["evt"] = masterArray["evt"][mask]

                for key,value in masterArray.items():
                    #if ( key in plottedVars or key == "finalweight" or key=="evt") \
                    #    and (len(mask)==len(value)):
                    #skimArray[key] = value[mask] #already applied
                    #which is correct below??
                    #if ( key in plottedVars or key == "finalweight" or key=="evt"):
                    if ( key in plottedVars or key == "finalweight"):
                        skimArray[key] = value

                print("skim after cuts ",len(skimArray["finalweight"]))


                #skimArrayPerCat[systematic+":"+cat+":"+processObj.nickname+":"+process] = skimArray
            #if process level ... end of cat loop
        #end of process loop?

            if process in ["FF_1","FF_2","FF_12"]: continue # rigging the FF multiple processes don't write out empty ones 
            filename = systematic+"_"+cat+"_"+processObj.nickname+"_"+process
            print("writing output file: ",filename)
            print("final yeild: " ,np.sum(masterArray["finalweight"]))
            try:
                os.mkdir("massOutputDir_"+args.outname)
            except:
                pass
            dataTypes =[[],[]]
            for branch in skimArray.keys():
                dataTypes[0].append(branch)
                dataTypes[1].append(skimArray[branch].dtype)


            randomkey = list(skimArray)[0]
            data = np.zeros(len(skimArray[randomkey]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
            fileout = ROOT.TFile.Open("massOutputDir_"+args.outname+"/"+filename,"recreate")
            fileout.cd()
            for branch in data.dtype.names:
                print(branch)
                try:
                    if len(skimArray[branch].shape) == 1:   # flat array important for jagged arrays of input data
                        data[branch] = skimArray[branch]
                    else:
                        data[branch] = skimArray[branch][:,0]
                except Exception as ex:
                    print(ex)
                    #continue

            root_numpy.array2tree(data,name=filename)
            fileout.Write()
            fileout.Close()
            del data
            del fileout
            del dataTypes
    #return skimArrayPerCat
    return


def checkRootFiles(systematics, allcats, processes,
                     finalDistributions, rootfiledir,
                     channel, outputstring):
   import os
   import glob
   import uproot
   import awkward as ak
   rootFiles = {}
   #print "the final distributions ",finalDistributions
   #cat = "mmmt_inclusive"
   finalSkims={}
   import shutil
   #list comprehension example:
   try:
       systematics[systematics.index("Events")]="Nominal"
   except:
       print("events tree not present")
   print("checking root files for sytematics ",systematics)
   finalSkims = {s:{c: dict() for c in allcats.keys()} for s in systematics}


   rootFiles={sys:glob.glob(rootfiledir+"/"+sys+"_*") for sys in systematics}


   mainTList = {}
   mainOutputTree={}
   emptyFiles={}
   nonemptyFiles={}
   #for sys, globfiles in rootFiles.items():
   #   for globfile in globfiles:
   #      #print "working on glob ",globfile,"  size   ",os.path.getsize(globfile) 
   #      for cat,catObj in allcats.items():
   #         for nickname, processObj in processes.items():
   #             for process in processObj.cuts.keys():
   #                 if str(globfile).split("/")[1]==sys+"_"+cat+"_"+nickname+"_"+process:
   #                     #print "found match ",sys+"_"+cat+"_"+nickname+"_"+process
   #                     fin = uproot.open(globfile)
   #                     if len(fin[sys+"_"+cat+"_"+nickname+"_"+process].keys())==12:
   #                         nonemptyFiles[sys+":"+cat+":"+nickname+":"+process] = processObj
   #                     else: 
   #                         emptyFiles[sys+":"+cat+":"+nickname+":"+process] = processObj
   resub=[]
   for sys in systematics:
       for cat,catObj in allcats.items():
           for nickname, processObj in processes.items():
               for process in processObj.cuts.keys():
                   for sys, globfiles in rootFiles.items():
                       found = 0
                       for globfile in globfiles:
                           if str(globfile).split("/")[1]==sys+"_"+cat+"_"+nickname+"_"+process:
                               found=1
                               fin = uproot.open(globfile)
                               #checking to see if there are output branches if so ... I assume it ran
                               if len(fin[sys+"_"+cat+"_"+nickname+"_"+process].keys())>2:
                                   nonemptyFiles[sys+":"+cat+":"+nickname+":"+process] = processObj
                               else: 
                                   emptyFiles[sys+":"+cat+":"+nickname+":"+process] = processObj
                                   resub.append(sys+":"+cat+":"+nickname+":"+process+"\n")
                       if not found:
                           emptyFiles[sys+":"+cat+":"+nickname+":"+process] = processObj
                           resub.append(sys+":"+cat+":"+nickname+":"+process+"\n")

                            
                        


   open("missingOutputFiles_"+args.channel+"_"+args.year+"_"+args.outname+".sh",'w').writelines(resub)
   return nonemptyFiles,emptyFiles

def slimskim(process,allcats,weightHistoDict,systematic):

    skimArrayPerSysCats={}
    #print "working on systematic ",systematic
    work_dict = {}
    #fin = uproot.open(process.file)
    with uproot.open(process.file) as fin:

        try:
            tree = fin[systematic]
        except:
            return skimArrayPerSysCats

        if systematic!="Events":
            syst_names = set(fin[systematic].keys())
            nom_names = set(fin["Events"].keys()) - syst_names
            work_dict.update(fin[systematic].arrays(list(syst_names)))
            work_dict.update(fin["Events"].arrays(list(nom_names)))
            skimArrayPerSysCats.update(makeCutsOnTreeArray(process,work_dict,allcats,weightHistoDict,systematic))
        else:
            skimArrayPerSysCats.update(makeCutsOnTreeArray(process,tree.arrays(library="np"),allcats,weightHistoDict,"Nominal"))

    del work_dict
    del tree
    return skimArrayPerSysCats




def slimskimoutput(process,allcats,weightHistoDict,systematic,massoutputdir,datadrivenPackage,jetWeightMultiplicity,fakefactorObj,args):

    skimArrayPerSysCats={}
    print("working on systematic ",systematic," opening file ",process.file)
    work_dict = {}
    with uproot.open(process.file) as fin:
        try:
           tree = fin[systematic]
        except:
           return

        if systematic!="Events":
            syst_names = set(fin[systematic].keys())
            nom_names = set(fin["Events"].keys()) - syst_names
            work_dict.update(fin[systematic].arrays(list(syst_names),library="np"))
            work_dict.update(fin["Events"].arrays(list(nom_names),library="np"))
            #skimArrayPerSysCats.update(makeCutsOnTreeArray(process,work_dict,allcats,weightHistoDict,systematic,jetWeightMultiplicity,fakefactorObj,args))
            makeCutsOnTreeArray(process,work_dict,allcats,weightHistoDict,systematic,jetWeightMultiplicity,fakefactorObj,args)
        else:
            #skimArrayPerSysCats.update(makeCutsOnTreeArray(process,tree.arrays(library="np"),allcats,weightHistoDict,"Nominal",jetWeightMultiplicity,fakefactorObj,args))
            makeCutsOnTreeArray(process,tree.arrays(library="np"),allcats,weightHistoDict,"Nominal",jetWeightMultiplicity,fakefactorObj,args)


        #createSlimOutput(skimArrayPerSysCats,massoutputdir)
    del skimArrayPerSysCats
    del work_dict
    del tree
    return

def createSlimOutput(skimArrayPerSysCats,outputdir):
   try:
       os.mkdir(outputdir)
   except:
       pass

       #print "dir exists " , outputdir
   print(skimArrayPerSysCats.keys())
   print("events ",len(skimArrayPerSysCats[skimArrayPerSysCats.keys()[0]].values())),
   for key,dictionary in skimArrayPerSysCats.items():
      dataTypes =[[],[]]
      #try:
      for branch in dictionary.keys():
          #if branch=="nickname":
              #print branch
              #print dictionary[branch].dtype
          dataTypes[0].append(branch)
          dataTypes[1].append(dictionary[branch].dtype)


      #data = np.zeros(len(dictionary[branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
      randomkey = dictionary.keys()[0]
      data = np.zeros(len(dictionary[randomkey]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
      #except:
      #    #data = np.zeros(0)
      #    continue

      filename = key.replace(":","_")
      fileout = ROOT.TFile.Open(outputdir+"/"+filename,"recreate")
      fileout.cd()
      for branch in data.dtype.names:
          if len(dictionary[branch].shape) == 1:   # flat array important for jagged arrays of input data
              data[branch] = dictionary[branch]
          else:
              data[branch] = dictionary[branch][:,0]

      #skimArrayPerCat[systematic+":"+cat+":"+process.nickname+":"+process] = skimArray
      #name = key.split(":")
      #root_numpy.array2tree(data,name=key.split(":")[0]+"_"+key.split(":")[-1])
      root_numpy.array2tree(data,name=filename)
      fileout.Write()
      fileout.Close()
      del data
      del fileout
   return




def combineRootFiles(systematics, allcats, processes,
                     finalDistributions, rootfiledir,
                     channel, outputstring):
   import os
   import glob
   import uproot
   import awkward as ak
   rootFiles = {}
   #print "the final distributions ",finalDistributions
   #cat = "mmmt_inclusive"
   finalSkims={}
   import shutil
   #list comprehension example:
   print("combining files for systematics ",systematics)
   try:
       systematics[systematics.index("Events")]="Nominal"
   except:
       print("events tree not present")
   finalSkims = {s:{c: dict() for c in allcats.keys()} for s in systematics}


   rootFiles={sys:glob.glob(rootfiledir+"/"+sys+"_*") for sys in systematics}


   mainTList = {}
   mainOutputTree={}
   for sys, globfiles in rootFiles.items():
      for globfile in globfiles:
         #print "working on glob ",globfile
         #process = globfile.split("_")[-1] # problem with data_obs or fake_W
         #process = globfile.split(sys)[0].split(channel+"_inclusive_")[1]
         #added loop for better processing?
         for cat,catObj in allcats.items():
             #if cat==args.channel+"_inclusive":
            for nickname, processObj in processes.items():
                for process in processObj.cuts.keys():
                    #if nickname in str(globfile) and nickname in str(globfile):
                    if str(globfile).split("/")[1]==sys+"_"+cat+"_"+nickname+"_"+process:
                        #print "found match ",sys+"_"+cat+"_"+nickname+"_"+process
                        with uproot.open(globfile) as fin:
                            tree = fin[sys+"_"+cat+"_"+nickname+"_"+process]
                            mainArrays = tree.arrays(library="np")
                            # mainArrays={}
                            # for branch,data in tree.items():
                            #     mainArrays[branch]=data.array(library="np")
                            #mainArrays["nickname"]=mainArrays["nickname"].array(library="np")
                            #print "tree ",sys+"_"+process," entries ",len(mainArrays["mll"])
                            for catDist, final in finalDistributions.items():
                                for processOut in final:
                                    if (processOut==process) and (catDist not in finalSkims[sys][cat]):
                                        #print "first output for process ",process," finalDist cat ",catDist
                                        finalSkims[sys][cat][catDist] = mainArrays
                                        continue
                                    elif (processOut==process) and (catDist in finalSkims[sys][cat]):
                                        #print "adding to finalskims ", catDist,"  for process ",process," finalDist cat ",catDist
                                        for branch in finalSkims[sys][cat][catDist].keys():
                                            if branch=="mask": continue
                                            if branch=="nickname":
                                                #print "changing nickname branch  "
                                                # print "nickname branch first ele ",finalSkims[sys][cat][catDist][branch][0]
                                                # print "nickname branch whole     ",finalSkims[sys][cat][catDist][branch]
                                                finalSkims[sys][cat][catDist][branch]=finalSkims[sys][cat][catDist][branch].astype("S40")
                                            #if finalSkims[sys][cat][catDist][branch].dtype == object:
                                            #    finalSkims[sys][cat][catDist][branch].astype('S')
                                            try:
                                                finalSkims[sys][cat][catDist][branch]=np.concatenate((finalSkims[sys][cat][catDist][branch],mainArrays[branch]))
                                            except:
                                                print("possible missing branch ",branch)
                                    else:
                                        continue

   #print "final skims example ", finalSkims["Nominal"]["mmtt_inclusive"]["rareBkg"]
   skimFile = ROOT.TFile("skimmed_"+outputstring+".root","recreate")
   skimFile.cd()
   for cat, catObj in allcats.items():
       #localtest="allsys"
       #skimFile = ROOT.TFile("skimmed_"+localtest+"_"+cat+".root","recreate")
       skimFile.cd()
       skimFile.mkdir(cat)
       skimFile.cd(cat)
       for sys in finalSkims.keys():
           #print "systematic ",sys
           dataTypes =[[],[]]
           #print finalSkims[sys][cat].values()
           #print "combining ",sys, " cat ",cat
           random_sample = list(finalSkims[sys][cat].values())[0]
           for branch in random_sample.keys():
               print(branch)
               if branch=="mask": continue
               #print random_sample[branch].dtype
               if branch=="nickname":
                   branchdatatype = "S40"
               else:
                   branchdatatype = random_sample[branch].dtype
               dataTypes[0].append(branch)
               dataTypes[1].append(branchdatatype)
           for catDist in finalSkims[sys][cat].keys():
               if branch=="mask": continue
               #print "on the final dist ",catDist
               data = np.zeros(len(finalSkims[sys][cat][catDist][branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
               for branch in data.dtype.names:
                   #print "working on branch ",branch
                   #print "branch datatype ",type(finalSkims[sys][cat][catDist][branch])
                   if len(finalSkims[sys][cat][catDist][branch].shape) == 1:   # flat array important for jagged arrays of input data
                       data[branch] = finalSkims[sys][cat][catDist][branch]
                   else:
                       data[branch] = finalSkims[sys][cat][catDist][branch][:,0]
               treeOut = root_numpy.array2tree(data, name=sys+"_"+catDist.split(";")[0])
               treeOut.Write()
   skimFile.Close()

   return 1


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


    print(args)
    allcats={}
    HAA_processes={}
    finalDistributions={}
    filelist = {}
    weightHistoDict={}
    jetWeightMultiplicity={}
    EventWeights={}
    datadrivenPackage={}


    allcats,HAA_processes,finalDistributions,\
    weightHistoDict,jetWeightMultiplicity,datadrivenPackage,fakefactorObj = initialize(args)


    #print(allcats)
    info('main line')
    #print(f)
    nums=[]
    skims=[]
    #skims={}
    payloadsdict={}
    payloads=[]

    if args.systematics:
        #depreciated
        #systematics =[ "scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
        #               "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
        #               "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
        #               "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"]
        #ultra legacy new ... all wps 
        systematics = [
        'scale_tUp',
        'scale_mUp',
        'scale_eUp',
        #mmtt
        'scale_t_1prong_TauLoose_MuoVLoose_EleVLooseUp',
        'scale_t_1prong1pizero_TauLoose_MuoVLoose_EleVLooseUp',
        'scale_t_3prong_TauLoose_MuoVLoose_EleVLooseUp',
        'scale_t_3prong1pizero_TauLoose_MuoVLoose_EleVLooseUp',
        #mmmt
        'scale_t_1prong_TauLoose_MuoTight_EleVLooseUp',
        'scale_t_1prong1pizero_TauLoose_MuoTight_EleVLooseUp',
        'scale_t_3prong_TauLoose_MuoTight_EleVLooseUp',
        'scale_t_3prong1pizero_TauLoose_MuoTight_EleVLooseUp',
        #mmet
        'scale_t_1prong_TauMedium_MuoVLoose_EleVTightUp',
        'scale_t_1prong1pizero_TauMedium_MuoVLoose_EleVTightUp',
        'scale_t_3prong_TauMedium_MuoVLoose_EleVTightUp',
        'scale_t_3prong1pizero_TauMedium_MuoVLoose_EleVTightUp',
        'scale_tDown',
        'scale_mDown',
        'scale_eDown',
        #mmtt
        'scale_t_1prong_TauLoose_MuoVLoose_EleVLooseDown',
        'scale_t_1prong1pizero_TauLoose_MuoVLoose_EleVLooseDown',
        'scale_t_3prong_TauLoose_MuoVLoose_EleVLooseDown',
        'scale_t_3prong1pizero_TauLoose_MuoVLoose_EleVLooseDown',
        #mmmt
        'scale_t_1prong_TauLoose_MuoTight_EleVLooseDown',
        'scale_t_1prong1pizero_TauLoose_MuoTight_EleVLooseDown',
        'scale_t_3prong_TauLoose_MuoTight_EleVLooseDown',
        'scale_t_3prong1pizero_TauLoose_MuoTight_EleVLooseDown',
        #mmet
        'scale_t_1prong_TauMedium_MuoVLoose_EleVTightDown',
        'scale_t_1prong1pizero_TauMedium_MuoVLoose_EleVTightDown',
        'scale_t_3prong_TauMedium_MuoVLoose_EleVTightDown',
        'scale_t_3prong1pizero_TauMedium_MuoVLoose_EleVTightDown',
        ]


        #takes too long to run FF closure region on systematics ... come on 
        allcats.pop(args.channel+"_FF_SS_validation")
    else:
        systematics =["Events"]
    #systematics =[ "Events"]

    jet_exc_samples = jet_exc_samples[args.year]
    nonemptyFiles,emptyFiles=\
    checkRootFiles(systematics, allcats, HAA_processes,\
                     finalDistributions, "massOutputDir_"+args.outname,\
                     args.channel, args.ffout)

    if args.rerun and not args.debug:
        #nonemptyFiles,emptyFiles=\
        #checkRootFiles(systematics, allcats, HAA_processes,\
        #                 finalDistributions, "massOutputDir_"+args.outname,\
        #                 args.channel, args.ffout)
        payloads=[]
        for syscatnickpro in emptyFiles.keys(): 
            sys = syscatnickpro.split(":")[0]
            if sys=="Nominal": sys="Events"
            payloads.append(
                (emptyFiles[syscatnickpro],allcats,
                weightHistoDict,sys,
                "massOutputDir_"+args.outname,
                 datadrivenPackage,jetWeightMultiplicity,fakefactorObj,args))

        print("MAKE SURE ALL CUTS ARE THE SAME")
        print("empty files ",emptyFiles.keys())
        #for payload in payloads:
        #    slimskimstar(payload)
        m = mp.Manager()
        logger_q = m.Queue()
        parallelable_data = [(1, logger_q), (2, logger_q)]
        pool  = mp.Pool(args.processingcores)

        pool.map(slimskimstar,payloads)
        pool.close()
        pool.join()
            
    if args.rerun and args.debug:
        payloads=[]
        for syscatnickpro in emptyFiles.keys(): 
            sys = syscatnickpro.split(":")[0]
            if sys=="Nominal": sys="Events"
            payloads.append(
                (emptyFiles[syscatnickpro],allcats,
                weightHistoDict,sys,
                "massOutputDir_"+args.outname,
                 datadrivenPackage,jetWeightMultiplicity,fakefactorObj,args))

        print("MAKE SURE ALL CUTS ARE THE SAME")
        print("empty files ",emptyFiles.keys())
        for payload in payloads:
            print("working on payload ", payload[0].nickname)
            slimskimstar(payload)
    if not args.rerun and not args.runprocesses:
        if not args.combine:
            for nickname, process in HAA_processes.items():
                for sys in systematics:
                    if sys=="Nominal": sys="Events"
                    payloads.append(
                        (process,allcats,
                        weightHistoDict,sys,
                        "massOutputDir_"+args.outname,
                         datadrivenPackage,jetWeightMultiplicity,fakefactorObj,args))
            #print(" PAYLOADS   ",payloads)

            if not args.debug:
                m = mp.Manager()
                logger_q = m.Queue()
                parallelable_data = [(1, logger_q), (2, logger_q)]
                pool  = mp.Pool(args.processingcores)

                pool.map(slimskimstar,payloads)
                pool.close()
                pool.join()
                #while not logger_q.empty():
                    #print logger_q.get()
            else:
                for payload in payloads:
                    print("working on payload ", payload[0].nickname)
                    #if "data" not in payload[0].nickname: continue
                    #if "FF" not in payload[0].nickname: continue
                    #print(payload[0].cuts)
                    #if payload[0].nickname!="SS_relaxed_data": continue
                    #if payload[0].nickname!="ZZTo4L": continue
                    slimskimstar(payload)
        


        print("root file generation computation time")
        print(datetime.datetime.now() - begin_time)
    if args.runprocesses:
        payloads=[]
        rps=[]
        for line in open(args.runprocesses,'r').readlines(): rps.append(line.strip())
        print("running on user input processes ",rps)
        for nickname, process in HAA_processes.items():
            if nickname not in rps: continue
            for sys in systematics:
                if sys=="Nominal": sys="Events"
                print("working on file ",process.file)
                print("working on distributions ",process.cuts.keys())
                payloads.append(
                    (process,allcats,
                    weightHistoDict,sys,
                    "massOutputDir_"+args.outname,
                    datadrivenPackage,jetWeightMultiplicity,fakefactorObj,args))
        for payload in payloads:
            print("working on payload ", payload[0].nickname)
            slimskimstar(payload)
            


    if args.combine:
        combineRootFiles(systematics, allcats, HAA_processes,
                     finalDistributions, "massOutputDir_"+args.outname,
                     args.channel, args.outname)
        print("combination successful")

        #shutil.rmtree("massOutputDir_"+args.outname)


    try:
        datadrivenPackage["fakemeasurefile"].Close()
    except:
        print("problem closing fake measurement helper file ... shutting down anyway")

    print("computation time")
    print(datetime.datetime.now() - begin_time)
    print("arguments used")
    print(args)
