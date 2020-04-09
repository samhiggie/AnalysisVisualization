#########################
#Author: Sam Higginbotham
'''

* File Name :Categories.py 

* Purpose : The Categories in the analysis... typically a set of extra cuts that define a phase space region of interest and stores the variables of interest in the fit and visualization

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
import ROOT
from Parametrization import Category
from functions import *

Inclusive = Category()
Inclusive.name = "inclusive"
Inclusive.variables = [["m_sv"]]
Inclusive.cuts = "mt_1<50.0"   #stuff like this
Inclusive.binning = [[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0]]


vbftest = Category()
vbftest.name = "vbf"
vbftest.variables = [["m_sv","mjj"],["njets"]]
vbftest.cuts = "mt_1>50.0&&njets>=2"   #stuff like this
vbftest.binning = [[[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0],[0.0,50.0,100.0,200.0,300.0,500.0]],[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]]


vbftestmela = Category()
vbftestmela.name = "vbfmela"
vbftestmela.variables = [["m_sv","mjj","mela"],["njets"]]
vbftestmela.cuts = "mt_1>50.0&&njets>=2"   #stuff like this
vbftestmela.binning = [[[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0],[0.0,50.0,100.0,200.0,300.0,500.0],[0.0,4.0,5.0]],[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0]]


#With HAA front end ... MUST select category -> channel ... 
#number = { 'eeet':1, 'eemt':2, 'eett':3, 'eeem':4, 'mmet':5, 'mmmt':6, 'mmtt':7, 'mmem':8, 'et':9, 'mt':10, 'tt':11 }
HAA_Inc_mmmt = Category()
HAA_Inc_mmmt.name = "mmmt_inclusive"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
#[["EQT"],["q_3","q_4"],"mult",">",0],[["EQT"],["q_1","q_2"],"mult",">",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #

HAA_Inc_mmmt.cuts["categoryCuts"]= [
["cat","==",6],["iso_3","<=",0.15]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],
#["mediumId_3",">=",1],["idDeepTau2017v2p1VSjet_4",">=",15.]]   #

#Trigger bit mapping... bits = [e.HLT_Ele27_eta2p1_WPTight_Gsf, e.HLT_Ele25_eta2p1_WPTight_Gsf, e.HLT_IsoMu24, e.HLT_IsoTkMu24, e.HLT_IsoMu27,e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,e.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ]
#2 -> 4 IsoMu24 ,  3 -> 8 IsoTkMu24  , 4 -> 16 IsoMu27  , 5 -> 32  e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ 
#HAA_Inc_mmmt.cuts["trigger"]= [[["triggerWord","band",4],"||",["triggerWord","band",8],"||",["triggerWord","band",16],"||",["triggerWord","band",32]]]   #
#if the statements are nested like this ... they are ors 
HAA_Inc_mmmt.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

#new variables that are defined on the fly. These are not in the Tree, but are functions of variables in the tree
#see function definitions above
HAA_Inc_mmmt.newvariables["mll-mtt"] = ["minus",["mll","m_vis"]]
HAA_Inc_mmmt.newvariablesbins = [[-40.0,-30.0,-20.0,-10.0,0.0,10.0,20.0,30.0,40.0,50.0,60.0]]
#HAA_Inc_mmmt.newvariables["mll-mtt"] = tauSFTool.getSFvsPT(e.pt_3,e.gen_match_3)["mll","mtt"]

HAA_Inc_mmmt.variables = [["m_vis"],
["AMass"],
#["mll"],
["mll"],
["pt_1"],]
#["pt_2"],
#["pt_3"],
#["pt_4"]]
HAA_Inc_mmmt.binning = [[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0],
#[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0],
#[0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.0,220.0,230.0,240.0,250.0,260.0,270.0,280.0,290.0,300.0,310.0,320.,330.,340.,350.,360.,370.,380.,390.,400.],
[0.0,20.0,40.0,60.0,80.0,100.0,120.,140.,160.,180.,200.,220.0,240.0,260.0,280.0,300.0,320.,340.,360.,380.,400.],
#[5.0,7.0,9.0,10.,12.,14.,16.0,18.,20.0,22.,24.0,26.,28.,30.0,35.0,40.0,50.0,60.0,70.0,80.0,90.,100.,120.,140.,160.,180.,200.],
[0.0,10.,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.,100.,110.,120.,130.,140.],
[0.0,25.0,50.,75.,100.,125.,150.,175.,200.]]

#HAA_Inc_mmmt.extraplots
allcats=[]
allcats.append(HAA_Inc_mmmt)

