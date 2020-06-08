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
import copy


#With HAA front end ... MUST select category -> channel ... 
#number = { 'eeet':1, 'eemt':2, 'eett':3, 'eeem':4, 'mmet':5, 'mmmt':6, 'mmtt':7, 'mmem':8, 'et':9, 'mt':10, 'tt':11 }
HAA_Inc_mmmt = Category()
HAA_Inc_mmmt.name = "mmmt_inclusive"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
#[["EQT"],["q_3","q_4"],"mult",">",0],[["EQT"],["q_1","q_2"],"mult",">",0],
#[["OR"],["q_3",">",0],["q_4",">",0][["EQT"],["q_1","q_2"],"mult",">",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #

HAA_Inc_mmmt.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],
["iso_3","<=",0.15],["mediumId_3",">=",1],
[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],
["idDeepTau2017v2p1VSjet_4",">=",2.]]
#,["idDeepTau2017v2p1VSjet_4",">=",15.] this is cutting most events! 

#Trigger bit mapping... bits = [e.HLT_Ele27_eta2p1_WPTight_Gsf, e.HLT_Ele25_eta2p1_WPTight_Gsf, e.HLT_IsoMu24, e.HLT_IsoTkMu24, e.HLT_IsoMu27,e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, e.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,e.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ]
#2 -> 4 IsoMu24 ,  3 -> 8 IsoTkMu24  , 4 -> 16 IsoMu27  , 5 -> 32  e.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ 
#HAA_Inc_mmmt.cuts["trigger"]= [[["triggerWord","band",4],"||",["triggerWord","band",8],"||",["triggerWord","band",16],"||",["triggerWord","band",32]]]   #
#if the statements are nested like this ... they are ors 
HAA_Inc_mmmt.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

#new variables that are defined on the fly. These are not in the Tree, but are functions of variables in the tree
#see function definitions above
HAA_Inc_mmmt.newvariables["mll-mtt"] = ["minus",[-40.0,-30.0,-20.0,-10.0,0.0,10.0,20.0,30.0,40.0,50.0,60.0],["mll","m_vis"],"[GeV]","mll - m_{vis}"]
HAA_Inc_mmmt.newvariables["charge_12"] = ["multi",[-2.0,-0.75,0.75,2.0],["q_1","q_2"],"","charge_1 #cross charge_2"]
HAA_Inc_mmmt.newvariables["charge_34"] = ["multi",[-2.0,-0.75,0.75,2.0],["q_3","q_4"],"","charge_3 #cross charge_4"]
#HAA_Inc_mmmt.newvariables["mll-mtt"] = tauSFTool.getSFvsPT(e.pt_3,e.gen_match_3)["mll","mtt"]

HAA_Inc_mmmt.vars = {
        #handle : ["root variable",[binning]or[[binning]],units,label]
        "m_vis":["m_vis",[50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0],"[GeV]","M_{vis}"],
        "AMass":["AMass",[100.0,120.,140.,160.,180.,200.,220.0,240.0,260.0,280.0,300.0,320.,340.,360.,380.,400.],"[GeV]","M_{4l Tot}"],
        "mll":["mll",[0.0,10.,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.,100.,110.,120.,130.,140.],"[GeV]","M_{ll}"],
        "pt_1":["pt_1",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#mu_{pt}"],
        "pt_2":["pt_2",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#mu_{pt}"],
        "pt_3":["pt_3",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#mu_{pt}"],
        "pt_4":["pt_4",[0.0,25.0,50.,75.,100.,125.,150.,175.,200.],"[GeV]","#tau_{pt}"],
        "pt_1_fine":["pt_1",[[40,0,200]],"[Gev]","P_{T}(#tau_{1})"],
        "eta_1":["eta_1",[[60,-3,3]],"","#eta(l_{1})"],
        "phi_1":["phi_1",[[60,-3,3]],"","#phi(l_{1})"],
        "iso_1":["iso_1",[[20,0,1]],"","relIso(l_{1})"],
        "dZ_1":["dZ_1",[[20,-0.2,0.2]],"[cm]","d_{z}(l_{1})"],
        "d0_1":["d0_1",[[20,-0.2,0.2]],"[cm]","d_{xy}(l_{1})"],
        "q_1":["q_1",[[3,-1.5,1.5]],"","charge(l_{1})"],
        "q_3":["q_3",[[3,-1.5,1.5]],"","charge(l_{3})"],
        "q_4":["q_4",[[3,-1.5,1.5]],"","charge(l_{4})"],
        "pt_2":["pt_2",[[40,0,200]],"[Gev]","P_{T}(l_{2})"],
        "eta_2":["eta_2",[[60,-3,3]],"","#eta(l_{2})"],
        "phi_2":["phi_2",[[60,-3,3]],"","#phi(l_{2})"],
        "iso_2":["iso_2",[[20,0,1]],"","relIso(l_{2})"],
        "dZ_2":["dZ_2",[[20,-0.2,0.2]],"[cm]","d_{z}(l_{2})"],
        "d0_2":["d0_2",[[20,-0.2,0.2]],"[cm]","d_{xy}(l_{2})"],
        "q_2":["q_2",[[3,-1.5,1.5]],"","charge(l_{2})"],
	    "iso_3":["iso_3",[[20,0,1]],"","relIso(l_{3})"],
        "pt_3":["pt_3",[[5,0,100]],"[Gev]","P_{T}(l_{3})"],
        "pt_3fine":["pt_3",[[40,0,100]],"[Gev]","P_{T}(l_{3})"],
        "eta_3":["eta_3",[[60,-3,3]],"","#eta(l_{3})"],
        "phi_3":["phi_3",[[60,-3,3]],"","#phi(l_{3})"],
        "dZ_3":["dZ_3",[[20,-0.2,0.2]],"[cm]","d_{z}(l_{3})"],
        "d0_3":["d0_3",[[20,-0.2,0.2]],"[cm]","d_{xy}(l_{3})"],
        "iso_4":["iso_4",[[20,0,1]],"","relIso(l_{4})"],
        "pt_4":["pt_4",[[5,0,100]],"[Gev]","P_{T}(l_{4})"],
        "pt_4fine":["pt_4",[[40,0,100]],"[Gev]","P_{T}(l_{4})"],
        "eta_4":["eta_4",[[60,-3,3]],"","#eta(l_{4})"],
        "phi_4":["phi_4",[[60,-3,3]],"","#phi(l_{4})"],
        "dZ_4":["dZ_4",[[20,-0.2,0.2]],"[cm]","d_{z}(l_{4})"],
        "d0_4":["d0_4",[[20,-0.2,0.2]],"[cm]","d_{xy}(l_{4})"],
        "njets":["njets",[[10,-0.5,9.5]],"","nJets"],
        "jpt_1":["jpt_1",[[60,0,300]],"[GeV]","Jet^{1} P_{T}"], 
        "jeta_1":["jeta_1",[[60,-3,3]],"","Jet^{1} #eta"],
        "jpt_2":["jpt_2",[[60,0,300]],"[GeV]","Jet^{2} P_{T}"], 
        "jeta_2":["jeta_2",[[6,-3,3]],"","Jet^{2} #eta"],
        "bpt_1":["bpt_1",[[40,0,200]],"[GeV]","BJet^{1} P_{T}"], 
        "bpt_2":["bpt_2",[[40,0,200]],"[GeV]","BJet^{2} P_{T}"], 
        "nbtag":["nbtag",[[5,-0.5,4.5]],"","nBTag"],
        "beta_1":["beta_1",[[60,-3,3]],"","BJet^{1} #eta"],
        "beta_2":["beta_2",[[60,-3,3]],"","BJet^{2} #eta"],
        "met":["met",[[20,0,200]],"[GeV]","#it{p}_{T}^{miss}"], 
        "met_phi":["met_phi",[[60,-3,3]],"","#it{p}_{T}^{miss} #phi"], 
        "puppi_phi":["puppi_phi",[[60,-3,3]],"","PUPPI#it{p}_{T}^{miss} #phi"], 
        "puppimet":["puppimet",[[20,0,200]],"[GeV]","#it{p}_{T}^{miss}"], 
        "mll_fine":["mll",[[40,20,100]],"[Gev]","m(l^{+}l^{-})"],
        "m_vis_fine":["m_vis",[[30,50,200]],"[Gev]","m(#tau#tau)"],
        "pt_tt":["pt_tt",[[40,0,200]],"[GeV]","P_{T}(#tau#tau)"],
        "A2_DR":["H_DR",[[60,0,6]],"","#Delta R(#tau#tau)"],
        "A2_tot":["H_tot",[[30,0,200]],"[GeV]","m_{T}tot(#tau#tau)"],
        "A1_Pt":["Z_Pt",[[60,0,300]],"[Gev]","P_T(l_{1}l_{2})"],
        "A1_DR":["Z_DR",[[60,0,6]],"[Gev]","#Delta R(l_{1}l_{2})"],
        "inTimeMuon_1":["inTimeMuon_1",[[3,-1.5,1.5]],"","inTimeMuon_1"],
        "isGlobal_1":["isGlobal_1",[[3,-1.5,1.5]],"","isGlobal_1"],
        "isTracker_1":["isTracker_1",[[3,-1.5,1.5]],"","isTracker_1"],
        "looseId_1":["looseId_1",[[3,-1.5,1.5]],"","looseId_1"],
        "mediumId_1":["mediumId_1",[[3,-1.5,1.5]],"","mediumId_1"],
        "Electron_mvaFall17V2noIso_WP90_1":["Electron_mvaFall17V2noIso_WP90_1",[[3,-1.5,1.5]],"","Electron_mvaFall17V2noIso_WP90_1"],
        "gen_match_1":["gen_match_1",[[30,-0.5,29.5]],"","gen_match_1"],
        "inTimeMuon_2":["inTimeMuon_2",[[3,-1.5,1.5]],"","inTimeMuon_2"],
        "isGlobal_2":["isGlobal_2",[[3,-1.5,1.5]],"","isGlobal_2"],
        "isTracker_2":["isTracker_2",[[3,-1.5,1.5]],"","isTracker_2"],
        "looseId_2":["looseId_2",[[3,-1.5,1.5]],"","looseId_2"],
        "mediumId_2":["mediumId_2",[[3,-1.5,1.5]],"","mediumId_2"],
        "Electron_mvaFall17V2noIso_WP90_2":["Electron_mvaFall17V2noIso_WP90_2",[[3,-1.5,1.5]],"","Electron_mvaFall17V2noIso_WP90_2"],
        "gen_match_2":["gen_match_2",[[30,-0.5,29.5]],"","gen_match_2"],
        "inTimeMuon_3":["inTimeMuon_3",[[3,-1.5,1.5]],"","inTimeMuon_3"],
        "isGlobal_3":["isGlobal_3",[[3,-1.5,1.5]],"","isGlobal_3"],
        "isTracker_3":["isTracker_3",[[3,-1.5,1.5]],"","isTracker_3"],
        "looseId_3":["looseId_3",[[3,-1.5,1.5]],"","looseId_3"],
        "mediumId_3":["mediumId_3",[[3,-1.5,1.5]],"","mediumId_3"],
        "Electron_mvaFall17V2noIso_WP90_3":["Electron_mvaFall17V2noIso_WP90_3",[[3,-1.5,1.5]],"","Electron_mvaFall17V2noIso_WP90_3"],
        "gen_match_3":["gen_match_3",[[30,-0.5,29.5]],"","gen_match_3"],
        "inTimeMuon_4":["inTimeMuon_4",[[3,-1.5,1.5]],"","inTimeMuon_4"],
        "isGlobal_4":["isGlobal_4",[[3,-1.5,1.5]],"","isGlobal_4"],
        "isTracker_4":["isTracker_4",[[3,-1.5,1.5]],"","isTracker_4"],
        "looseId_4":["looseId_4",[[3,-1.5,1.5]],"","looseId_4"],
        "mediumId_4":["mediumId_4",[[3,-1.5,1.5]],"","mediumId_4"],
        "Electron_mvaFall17V2noIso_WP90_4":["Electron_mvaFall17V2noIso_WP90_4",[[3,-1.5,1.5]],"","Electron_mvaFall17V2noIso_WP90_4"],
        "gen_match_4":["gen_match_4",[[30,-0.5,29.5]],"","gen_match_4"],
}

HAA_Inc_mmmt_noIds = copy.deepcopy(HAA_Inc_mmmt)

HAA_Inc_mmmt_noIds.name = "mmmt_inclusive_noIds"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_noIds.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3]]

HAA_Inc_mmmt_noIds.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.]]


HAA_Inc_mmmt_charge = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_charge.name = "mmmt_inclusive_charge"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_charge.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0]]
HAA_Inc_mmmt_charge.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.]]


HAA_Inc_mmmt_iso = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_iso.name = "mmmt_inclusive_iso"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_iso.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_iso.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15]]
#["mediumId_3",">=",1],
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#,["idDeepTau2017v2p1VSjet_4",">=",15.] this is cutting most events! 
HAA_Inc_mmmt_iso.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #



HAA_Inc_mmmt_muon = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_muon.name = "mmmt_inclusive_muon"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_muon.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_muon.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],
[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#,["idDeepTau2017v2p1VSjet_4",">=",15.] this is cutting most events! 
HAA_Inc_mmmt_muon.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

HAA_Inc_mmmt_tauid = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_tauid.name = "mmmt_inclusive_tauid"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_tauid.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_tauid.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],["idDecayModeNewDMs_4",">",0.]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",15.]this is cutting most events! 
HAA_Inc_mmmt_tauid.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

HAA_Inc_mmmt_deeptauid_medium = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_deeptauid_medium.name = "mmmt_inclusive_deeptauid_medium"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_deeptauid_medium.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_deeptauid_medium.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],["idDeepTau2017v2p1VSjet_4",">=",16.]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",15.]this is cutting most events! 
HAA_Inc_mmmt_deeptauid_medium.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

HAA_Inc_mmmt_deeptauid_loose = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_deeptauid_loose.name = "mmmt_inclusive_deeptauid_loose"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_deeptauid_loose.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_deeptauid_loose.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],["idDeepTau2017v2p1VSjet_4",">=",8.]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",15.]this is cutting most events! 
HAA_Inc_mmmt_deeptauid_loose.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #
allcats=[]

HAA_Inc_mmmt_deeptauid_vloose = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_deeptauid_vloose.name = "mmmt_inclusive_deeptauid_vloose"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_deeptauid_vloose.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_deeptauid_vloose.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],["idDeepTau2017v2p1VSjet_4",">=",4.]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",15.]this is cutting most events! 
HAA_Inc_mmmt_deeptauid_vloose.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #

HAA_Inc_mmmt_deeptauid_vvloose = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_deeptauid_vvloose.name = "mmmt_inclusive_deeptauid_vvloose"
#all things in brackets are "annnd"s a single cut can then be "OR" also even an equation "EQT" like q_1*q_2
HAA_Inc_mmmt_deeptauid_vvloose.cuts["preselection"]= [
#["isTrig_1","!=",0],
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_deeptauid_vvloose.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],["iso_3","<=",0.15],
["mediumId_3",">=",1],["idDeepTau2017v2p1VSjet_4",">=",2.]]
#[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",15.]this is cutting most events! 
HAA_Inc_mmmt_deeptauid_vvloose.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]
#[["IF"],[["isTrig_1","==",1],["isTrig_2","==",0]],["THEN"],[["pt_1",">=",28]]]]   #




#override for now skip plotting other variables
#HAA_Inc_mmmt.vars = {
#        #handle : ["root variable",[binning]or[[binning]],units,label]
#        "AMass":["AMass",[100.0,120.,140.,160.,180.,200.,220.0,240.0,260.0,280.0,300.0,320.,340.,360.,380.,400.],"[GeV]","M_{4l Tot}"],
#}

#for data driven method... 
# adding new mt cut variables 
#HAA_Inc_mmmt.newvariables["mll-mtt"] = ["minus",[-40.0,-30.0,-20.0,-10.0,0.0,10.0,20.0,30.0,40.0,50.0,60.0],["mll","m_vis"],"[GeV]","mll - m_{vis}"]
#HAA_Inc_mmmt.newvariables["ptMiss"] = ["TLorentz",[-40.0,-30.0,-20.0,-10.0,0.0,10.0,20.0,30.0,40.0,50.0,60.0],["pt_3","eta_3","phi_3",0.000511],"[GeV]","mll - m_{vis}","multiob"]
#HAA_Inc_mmmt.newvariables["mut"] = ["TLorentz",[-40.0,-30.0,-20.0,-10.0,0.0,10.0,20.0,30.0,40.0,50.0,60.0],["pt_3","eta_3","phi_3",0.000511],"[GeV]","mll - m_{vis}","multiob"]
#
#HAA_Inc_mmmt.newvariables["muMET"] = ["addTLorentz",[-40.0,-30.0,-20.0,-10.0,0.0,10.0,20.0,30.0,40.0,50.0,60.0],["mut","ptMiss"],"[GeV]","","multiob"]
#
#HAA_Inc_mmmt.newvariables["muMETmt"] = ["TLorentzMT",[-40.0,-30.0,-20.0,-10.0,0.0,10.0,20.0,30.0,40.0,50.0,60.0],["muMET"],"[GeV]",""]
#HAA_Inc_mmmt_FF_loose_SS.newvarcuts["mtcuts"]= [["muMETmt","<","40"]]

#HAA_Inc_mmmt.vars = {
#        "pt_3":["pt_3",[[5,0,100]],"[Gev]","P_{T}(l_{3})"],
#        "pt_4":["pt_4",[[5,0,100]],"[Gev]","P_{T}(l_{4})"],
#}

HAA_Inc_mmmt_FF_SS_1_tight = copy.deepcopy(HAA_Inc_mmmt)
HAA_Inc_mmmt_FF_SS_1_tight.name = "mmmt_inclusive_FF_SS_1_tight"
HAA_Inc_mmmt_FF_SS_1_tight.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult",">",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_FF_SS_1_tight.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],
["iso_3","<=",0.15],["mediumId_3",">=",1],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],
["idDeepTau2017v2p1VSjet_4",">=",2.]]
HAA_Inc_mmmt_FF_SS_1_tight.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]

#invert muon criteria
HAA_Inc_mmmt_FF_SS_1_loose = copy.deepcopy(HAA_Inc_mmmt_FF_SS_1_tight)
HAA_Inc_mmmt_FF_SS_1_loose.name = "mmmt_inclusive_FF_SS_1_loose"
HAA_Inc_mmmt_FF_SS_1_loose.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult",">",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_FF_SS_1_loose.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],
["iso_3",">",0.15],["mediumId_3","<",1],["isGlobal_3","<",1],["isTracker_3","<",1],
#["iso_3","<=",0.15],["mediumId_3",">=",1],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]
["idDeepTau2017v2p1VSjet_4",">=",2.]]
HAA_Inc_mmmt_FF_SS_1_loose.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]

HAA_Inc_mmmt_FF_SS_2_tight = copy.deepcopy(HAA_Inc_mmmt_FF_SS_1_tight)
HAA_Inc_mmmt_FF_SS_2_tight.name = "mmmt_inclusive_FF_SS_2_tight"
HAA_Inc_mmmt_FF_SS_2_tight.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult",">",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_FF_SS_2_tight.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],
["iso_3","<=",0.15],["mediumId_3",">=",1],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],
["idDeepTau2017v2p1VSjet_4",">=",2.]]
HAA_Inc_mmmt_FF_SS_2_tight.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]

#invert tau criteria
HAA_Inc_mmmt_FF_SS_2_loose = copy.deepcopy(HAA_Inc_mmmt_FF_SS_1_tight)
HAA_Inc_mmmt_FF_SS_2_loose.name = "mmmt_inclusive_FF_SS_2_loose"
HAA_Inc_mmmt_FF_SS_2_loose.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult",">",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_FF_SS_2_loose.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],
["iso_3","<=",0.15],["mediumId_3",">=",1],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]],
["idDeepTau2017v2p1VSjet_4","<",2.]]
HAA_Inc_mmmt_FF_SS_2_loose.cuts["trigger"]=[
[["OR"],["triggerWord","band",4],["triggerWord","band",8],["triggerWord","band",16],["triggerWord","band",32]]]


#HAA_Inc_mmmt_FF_loose_SS = copy.deepcopy(HAA_Inc_mmmt_FF_tight_SS)
#HAA_Inc_mmmt_FF_loose_SS.cuts["preselection"]= [
#["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
#["eta_4","absl",2.3],
#[["EQT"],["q_3","q_4"],"mult",">",0],[["EQT"],["q_1","q_2"],"mult","<",0],
#["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
#HAA_Inc_mmmt_FF_loose_SS.cuts["categoryCuts"]= [
#["cat","==",6],["AMass",">=",120.]]

HAA_Inc_mmmt_FF_OS_1 = copy.deepcopy(HAA_Inc_mmmt_FF_SS_1_tight)
HAA_Inc_mmmt_FF_OS_1.name = "mmmt_inclusive_FF_OS_1"
HAA_Inc_mmmt_FF_OS_1.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_FF_OS_1.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],
#["iso_3","<=",0.15],["mediumId_3",">=",1],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]
["idDeepTau2017v2p1VSjet_4",">=",2.]]

HAA_Inc_mmmt_FF_OS_2 = copy.deepcopy(HAA_Inc_mmmt_FF_OS_1)
HAA_Inc_mmmt_FF_OS_2.name = "mmmt_inclusive_FF_OS_2"
HAA_Inc_mmmt_FF_OS_2.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_FF_OS_2.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.],
["iso_3","<=",0.15],["mediumId_3",">=",1],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",2.]]

HAA_Inc_mmmt_FF_OS_12 = copy.deepcopy(HAA_Inc_mmmt_FF_OS_1)
HAA_Inc_mmmt_FF_OS_12.name = "mmmt_inclusive_FF_OS_12"
HAA_Inc_mmmt_FF_OS_12.cuts["preselection"]= [
["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],
["eta_4","absl",2.3],
[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],
["iso_1","<=",0.2],["iso_2","<=",0.2]]   #
HAA_Inc_mmmt_FF_OS_12.cuts["categoryCuts"]= [
["cat","==",6],["AMass",">=",120.]]
#["iso_3","<=",0.15],["mediumId_3",">=",1],[["OR"],["isGlobal_3",">=",1],["isTracker_3",">=",1]]]
#["idDeepTau2017v2p1VSjet_4",">=",2.]]

ffbasecuts = [["pt_1",">",5.],["pt_2",">",5.],["pt_3",">",5.],["pt_4",">",18.5],["eta_4","absl",2.3],[["EQT"],["q_3","q_4"],"mult","<",0],[["EQT"],["q_1","q_2"],"mult","<",0],["iso_1","<=",0.2],["iso_2","<=",0.2],["cat","==",6],["AMass",">=",120.]]

#allcats.append(HAA_Inc_mmmt)
#current ff here may 8th 2020
#allcats.append(HAA_Inc_mmmt_FF_SS_1_tight)
#allcats.append(HAA_Inc_mmmt_FF_SS_1_loose)
#allcats.append(HAA_Inc_mmmt_FF_SS_2_tight)
#allcats.append(HAA_Inc_mmmt_FF_SS_2_loose)

allcats.append(HAA_Inc_mmmt_deeptauid_medium)
allcats.append(HAA_Inc_mmmt_deeptauid_loose)
allcats.append(HAA_Inc_mmmt_deeptauid_vloose)
allcats.append(HAA_Inc_mmmt_deeptauid_vvloose)


#allcats.append(HAA_Inc_mmmt_FF_SS_12)
#allcats.append(HAA_Inc_mmmt_FF_OS_1)
#allcats.append(HAA_Inc_mmmt_FF_OS_2)
#allcats.append(HAA_Inc_mmmt_FF_OS_12)

#allcats.append(HAA_Inc_mmmt_FF_loose_SS)
#allcats.append(HAA_Inc_mmmt_FF_tight_SS)
#allcats.append(HAA_Inc_mmmt_FF_loose_OS)
#allcats.append(HAA_Inc_mmmt_FF_tight_OS)

#allcats.append(HAA_Inc_mmmt_noIds)
#allcats.append(HAA_Inc_mmmt_charge)
#allcats.append(HAA_Inc_mmmt_iso)
#allcats.append(HAA_Inc_mmmt_muon)
#allcats.append(HAA_Inc_mmmt_deeptauid)
