#########################
#Author: Sam Higginbotham
'''

* File Name :

* Purpose :

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
#HTT
ProCuts = {}
ProCuts["TTT"]="gen_match_2>=5"
ProCuts["TTL"]="gen_match_2==5"
ProCuts["ZTT"]="gen_match_2>=5"
ProCuts["ZTL"]="gen_match_2==5"
ProCuts["ZL"]="gen_match_2<5"

ProCuts["jetFakes"]="HLT_Any>0"

HAA_ProCuts = {}

HAA_ProCuts["a15"]="16*00005."
HAA_ProCuts["a20"]="16*00005."
HAA_ProCuts["a25"]="16*00005."
HAA_ProCuts["a30"]="16*00005."
HAA_ProCuts["a35"]="16*00005."
HAA_ProCuts["a40"]="16*00005."
HAA_ProCuts["a45"]="16*00005."
HAA_ProCuts["a50"]="16*00005."
HAA_ProCuts["a55"]="16*00005."
HAA_ProCuts["a60"]="16*00005."
#HAA_ProCuts["TTT"]="gen_match_4>=5"
#HAA_ProCuts["ZZ"]="gen_match_4==5"
HAA_ProCuts["WZJ"]="(gen_match_3>5||gen_match_4>5)"
HAA_ProCuts["WJ"]="(gen_match_3>5||gen_match_4>5)"
HAA_ProCuts["ZJ"]="(gen_match_3>5||gen_match_4>5)"
HAA_ProCuts["WZL"]="(gen_match_3>=5||gen_match_4>=5)"
HAA_ProCuts["WL"]="(gen_match_3>=5||gen_match_4>=5)"
HAA_ProCuts["ZL"]="(gen_match_3>=5||gen_match_4>=5)"
HAA_ProCuts["ZTT"]="gen_match_4==5"
HAA_ProCuts["TTT"]="gen_match_4==5"
HAA_ProCuts["TTJ"]="gen_match_4>5"
HAA_ProCuts["TTL"]="gen_match_4>=5"
#HAA_ProCuts["EWK"]=""
