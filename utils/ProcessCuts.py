#########################
#Author: Sam Higginbotham
'''

* File Name :

* Purpose :

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
ProCuts = {}
ProCuts["preselection"]="pt_1>21&&pt_2>30&&npv>0&&id_m_medium_1>0&&iso_1<0.15&&byTightIsolationMVArun2v1DBoldDMwLT_2>0.5&&tightMuons<=1&&tightElectrons==0&&diLeptons==0&&againstMuonTight3_2>0&&againstElectronVLooseMVA6_2>0"
ProCuts["trigger"]="HLT_Any>0"
ProCuts["TTT"]="gen_match_2>=5"
ProCuts["TTL"]="gen_match_2==5"
ProCuts["ZTT"]="gen_match_2>=5"
ProCuts["ZTL"]="gen_match_2==5"
ProCuts["ZL"]="gen_match_2<5"

ProCuts["jetFakes"]="HLT_Any>0"
