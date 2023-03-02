#########################
#Author: Sam Higginbotham
'''

* File Name :Weights.py

* Purpose :Contains the mapping between weights and processes

* Creation Date : 05-02-2020

* Last Modified :

'''
#########################
CommonWeights={}
#CommonWeights["lumi2018"]=[59740]
CommonWeights["lumi2018"]=[59830]
#CommonWeights["lumi2016"]=[35900]#pb^-1
CommonWeights["lumi2016pre"]=[36310]#pb^-1
CommonWeights["lumi2016post"]=[36310]#pb^-1
#CommonWeights["lumi2017"]=[41800]#pb^-1
CommonWeights["lumi2017"]=[41480]#pb^-1

jet_exc_samples = {}

jet_exc_samples["2018"] = [
                        "DYJetsToLL",
                        "DYJetsToLLM10to50_ext1",
                        "DYJetsToLLM10to50",
                        "DY1JetsToLL" ,
                        "DY2JetsToLL" ,
                        "DY3JetsToLL" ,
                        "DY4JetsToLL" ,
                        "WJetsToLNu" ,
                        "W1JetsToLNu" ,
                        "W2JetsToLNu" ,
                        "W3JetsToLNu",
                        "W4JetsToLNu",
                        ]

jet_exc_samples["2017"] = [
                        "DYJetsToLL_ext1",
                        "DYJetsToLL",
                        "DYJetsToLLM10to50_ext1",
                        "DYJetsToLLM10to50",
                        "DY1JetsToLL" ,
                        "DY1JetsToLL_ext1" ,
                        "DY2JetsToLL" ,
                        "DY2JetsToLL_ext1" ,
                        "DY3JetsToLL" ,
                        "DY3JetsToLL_ext1" ,
                        "DY4JetsToLL" ,
                        "WJetsToLNu_ext1" ,
                        "W1JetsToLNu" ,
                        "W2JetsToLNu",
                        "W3JetsToLNu",
                        "W4JetsToLNu"
                        ]

jet_exc_samples["2016pre"] = ["DYJetsToLLext1",
                        "DYJetsToLLext2" ,
                        "DYJetsToLLM10to50" ,
                        "DY1JetsToLL" ,
                        "DY2JetsToLL" ,
                        "DY3JetsToLL" ,
                        "DY4JetsToLL" ,
                        "WJetsToLNu" ,
                        "WJetsToLNu_ext2" ,
                        "W1JetsToLNu" ,
                        "W2JetsToLNu_ext1" ,
                        "W3JetsToLNu",
                        "W3JetsToLNu_ext1",
                        "W4JetsToLNu",
                        "W4JetsToLNu_ext1",
                        "W4JetsToLNu_ext2"
                        ]
jet_exc_samples["2016post"] = ["DYJetsToLLext1",
                        "DYJetsToLLext2" ,
                        "DYJetsToLLM10to50" ,
                        "DY1JetsToLL" ,
                        "DY2JetsToLL" ,
                        "DY3JetsToLL" ,
                        "DY4JetsToLL" ,
                        "WJetsToLNu" ,
                        "WJetsToLNu_ext2" ,
                        "W1JetsToLNu" ,
                        "W2JetsToLNu_ext1" ,
                        "W3JetsToLNu",
                        "W3JetsToLNu_ext1",
                        "W4JetsToLNu",
                        "W4JetsToLNu_ext1",
                        "W4JetsToLNu_ext2"
                        ]

jetIncOnly = {}
jetIncOnly["2018"] = [["DYJetsToLL"],["WJetsToLNu"]]
jetIncOnly["2017"] = [["DYJetsToLL","DYJetsToLL_ext1"],["WJetsToLNu","WJetsToLNu_ext1"]]
jetIncOnly["2016pre"] = [["DYJetsToLLext1","DYJetsToLLext2"],["WJetsToLNu","WJetsToLNu_ext2"]]
jetIncOnly["2016post"] = [["DYJetsToLLext1","DYJetsToLLext2"],["WJetsToLNu","WJetsToLNu_ext2"]]


