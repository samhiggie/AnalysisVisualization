import ROOT
import numpy as np 
import uproot

def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

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
    from utils.Parametrization import fakefactor

    #Gather the analysis datasets and info
    sampleDict = {}
    csvfile = args.csvfile
    categories = args.categories
    processes = args.processes

    dir = args.dir



    for line in open(csvfile,'r').readlines() :
            # structure for csv to sampleDict conversion
            #[nickname]        = [category,xsec,numberOfEvents,finishedEvents,idk?,DASDataset]
            print "reading line "
            print line
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
        tempcat = Category()
        tempcat.name=categories[category]['name']
        tempcat.cuts=categories[category]['cuts']
        tempcat.newvariables=categories[category]['newvariables']
        tempcat.vars=categories[category]['vars']
        allcats[tempcat.name]=tempcat
    #loading standard processes
    HAA_processes={}

    nevents = 250000

    if not (args.datadrivenZH or args.datadrivenSM):
        for sample in sampleDict.keys():
            temppro = Process()
            temppro.nickname=sample
            temppro.file=dir+sample+"_"+args.year+".root"

            frooin = ROOT.TFile.Open(temppro.file,"read")
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue","genweight":"Generator_weight"}
            temppro.cuts={sampleDict[sample][0]:""}
            if "ggTo2mu2tau" in sample:
                if args.year==2016: nevents = 250000
                if args.year==2017: nevents = 350000
                if args.year==2017: nevents = 500000
                if args.extract:
                    temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(48.37*0.001)} # SM Higgs xsec [pb] x BR Haa x 5 for DataMC control plots

                else:
                    temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(137.5*31.05*0.00005)} # worked before
            HAA_processes[temppro.nickname]=temppro

    if (args.datadrivenZH or args.datadrivenSM):
        for sample in sampleDict.keys():
            temppro = Process()
            temppro.nickname=sample
            temppro.file=dir+sample+"_"+args.year+".root"
            temppro.weights={"xsec":sampleDict[sample][1],"nevents":sampleDict[sample][3],"PU":"weightPUtrue","genweight":"Generator_weight"}
            if args.channel=="mmtt":
                truetau = [
                            [["OR"],
                            ["gen_match_3","==",5],
                            ["gen_match_4","==",5]
                            ]]
            if args.channel=="mmem":
                truetau = [ [["OR"],
                            ["gen_match_3","==",15]
                            ],
                            [["OR"],
                            ["gen_match_4","==",15]
                            ] ]
            if args.channel=="mmmt" or args.channel=="mmet":
                truetau = [[["OR"],
                            ["gen_match_3","==",15],
                            ["gen_match_4","==",5]
                            ]]
            temppro.cuts={sampleDict[sample][0]:truetau} #ONLY SELECT PRMOPT FOR MC!
            if "ggTo2mu2tau" in sample:
                if args.year==2016: nevents = 250000
                if args.year==2017: nevents = 350000
                if args.year==2017: nevents = 500000
                #for visualization!
                if args.extract:
                    #temppro.weights={"xsec":1,"nevents":250000,"theoryXsec":(48.37*0.0001)} # SM Higgs xsec x BR Haa x 5 for DataMC control plots
                    temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(48.37*0.001)} # SM Higgs xsec x BR Haa x 5 for DataMC control plots

                else:
                    temppro.weights={"xsec":1,"nevents":nevents,"theoryXsec":(137.5*31.05*0.00005)} # worked before
            HAA_processes[temppro.nickname]=temppro

    if (args.datameasureZH):
        try:
            os.mkdir("FFhistos_"+str(args.ffout))
        except:
            print("directory exists")
        for proObj in HAA_processes.keys():
            if proObj!="data_obs" and proObj!="FF" and "ggTo2mu2tau" not in proObj:
                if args.channel=="mmmt" or args.channel=="mmet":
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
    elif (args.datadrivenZH or args.datadrivenSM):
        for proObj in HAA_processes.keys():
            if proObj!="data_obs" and proObj!="FF" and not args.skim:
                HAA_processes[proObj].cuts["fake1_"+str(proObj)] = [["gen_match_3","==",0]]
                HAA_processes[proObj].cuts["fake2_"+str(proObj)] = [["gen_match_4","==",0]]
                if "ggTo2mu2tau" in proObj:
                    print "skipping signal prompt mc"
                    continue  #this hasn't been added yet ... this is needed to omit signal!!
                HAA_processes[proObj].cuts["prompt1"] = [["gen_match_3","!=",0]]
                HAA_processes[proObj].cuts["prompt2"] = [["gen_match_4","!=",0]]


    #loading special processes ... fake factor and data!
    for process in processes_special:
        temppro = Process()
        temppro.nickname=processes_special[process]['nickname']
        temppro.cuts=processes_special[process]['cuts']
        temppro.weights=processes_special[process]['weights']
        temppro.file=dir+processes_special[process]['file']
        HAA_processes[temppro.nickname]=temppro
    #print "==========================================================================================="
    #print " All the processes !!!!!  ",HAA_processes
    #print "==========================================================================================="


    #precise weights for jet multiplicity
    from utils.Weights import jet_exc_samples
    from utils.Weights import jetIncOnly
    inc_samps = jetIncOnly[args.year]
    jet_exc_samples = jet_exc_samples[args.year]

    #file list for easy reading
    for proObj in HAA_processes.keys():
        filelist[proObj]=HAA_processes[proObj].file

    # adding kfactor to relevent samples
    for sample in sampleDict.keys():
        if sample in jet_exc_samples and sample.startswith("DY"):
            HAA_processes[sample].weights.update({"kfactor":1.1637})
        if sample in jet_exc_samples and sample.startswith("W"):
            HAA_processes[sample].weights.update({"kfactor":1.221})

    weightHistoDict = {}
    for nickname, filename in filelist.iteritems():
        #frooin = ROOT.TFile(dir+filename)
        frooin = ROOT.TFile(filename)
        hists = None
        try:
            if nickname.startswith("DYJetsToLL"):
                hists=[frooin.Get("hWeights").Clone(),
                       frooin.Get("DY1genWeights").Clone(),
                       frooin.Get("DY2genWeights").Clone(),
                       frooin.Get("DY3genWeights").Clone(),
                       frooin.Get("DY4genWeights").Clone(),
                       ]
            elif nickname.startswith("WJetsToLNu"):
                hists=[frooin.Get("hWeights").Clone(),
                      frooin.Get("W1genWeights").Clone(),
                      frooin.Get("W2genWeights").Clone(),
                      frooin.Get("W3genWeights").Clone(),
                      frooin.Get("W4genWeights").Clone(),
                      ]
            else:
                hists = frooin.Get("hWeights").Clone()
                hists.SetDirectory(0)
            if isinstance(hists, list):
                for hist in hists:
                    hist.SetDirectory(0)
            weightHistoDict[nickname] = hists
            frooin.Close()
        except:
            print nickname," hWeights prob doesn't exist?"

    print weightHistoDict
    #for key in weightHistoDict.keys():
        #if weightHistoDict[key].Get
    #print weightHistoDict["DYJetsToLLext1"]
    #print weightHistoDict["WJetsToLNu"]


    jetWeightMultiplicity = {}
    NJetWeights = {}
    # for all the DYJet and WJet files make sure that the exts are summed PRIOR to makeCutsOnTreeArray
    #print inc_samps
    for inc_group in inc_samps:
        sumOfWeights = np.zeros(5)
        for sample in inc_group:
            weights = weightHistoDict[sample]
            sumOfWeights += [w.GetSumOfWeights() for w in weights]
        sumOfWeights[sumOfWeights==0] = 1e-10 # if value is 0 to stop infs
        #print(sumOfWeights)
        for sample in inc_group:
            jetWeightMultiplicity[sample] = sumOfWeights


    if args.year=="2018":
        DY1JetsFile = ROOT.TFile.Open(filelist["DY1JetsToLL"],"read")
        jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
        DY2JetsFile = ROOT.TFile.Open(filelist["DY2JetsToLL"],"read")
        jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
        DY3JetsFile = ROOT.TFile.Open(filelist["DY3JetsToLL"],"read")
        jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
        DY4JetsFile = ROOT.TFile.Open(filelist["DY4JetsToLL"],"read")
        jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()

        W1JetsFile = ROOT.TFile.Open(filelist["W1JetsToLNu"],"read")
        jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()

        W2JetsFile = ROOT.TFile.Open(filelist["W2JetsToLNu"],"read")
        jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()

        W3JetsFile = ROOT.TFile.Open(filelist["W3JetsToLNu"],"read")
        jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()

        W4JetsFile = ROOT.TFile.Open(filelist["W4JetsToLNu"],"read")
        jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()

    if args.year=="2017":
        DY1JetsFile = ROOT.TFile.Open(filelist["DY1JetsToLL"],"read")
        jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
        DY1JetsFile_ext1 = ROOT.TFile.Open(filelist["DY1JetsToLL_ext1"],"read")
        jetWeightMultiplicity["DY1JetsToLL_ext1"]=DY1JetsFile_ext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["DY1JetsToLL_ext1"]+=jetWeightMultiplicity["DY1JetsToLL"]
        jetWeightMultiplicity["DY1JetsToLL"]+=jetWeightMultiplicity["DY1JetsToLL_ext1"]

        DY2JetsFile = ROOT.TFile.Open(filelist["DY2JetsToLL"],"read")
        jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
        DY2JetsFile_ext1 = ROOT.TFile.Open(filelist["DY2JetsToLL_ext1"],"read")
        jetWeightMultiplicity["DY2JetsToLL_ext1"]=DY2JetsFile_ext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["DY2JetsToLL_ext1"]+=jetWeightMultiplicity["DY2JetsToLL"]
        jetWeightMultiplicity["DY2JetsToLL"]+=jetWeightMultiplicity["DY2JetsToLL_ext1"]

        DY3JetsFile = ROOT.TFile.Open(filelist["DY3JetsToLL"],"read")
        jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
        DY3JetsFile_ext1 = ROOT.TFile.Open(filelist["DY3JetsToLL_ext1"],"read")
        jetWeightMultiplicity["DY3JetsToLL_ext1"]=DY3JetsFile_ext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["DY3JetsToLL_ext1"]+=jetWeightMultiplicity["DY3JetsToLL"]
        jetWeightMultiplicity["DY3JetsToLL"]+=jetWeightMultiplicity["DY3JetsToLL_ext1"]

        DY4JetsFile = ROOT.TFile.Open(filelist["DY4JetsToLL"],"read")
        jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()

        W1JetsFile = ROOT.TFile.Open(filelist["W1JetsToLNu"],"read")
        jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()

        W2JetsFile = ROOT.TFile.Open(filelist["W2JetsToLNu"],"read")
        jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()

        W3JetsFile = ROOT.TFile.Open(filelist["W3JetsToLNu"],"read")
        jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()

        W4JetsFile = ROOT.TFile.Open(filelist["W4JetsToLNu"],"read")
        jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()

    if args.year=="2016":
        DY1JetsFile = ROOT.TFile.Open(filelist["DY1JetsToLL"],"read")
        jetWeightMultiplicity["DY1JetsToLL"]=DY1JetsFile.Get("hWeights").GetSumOfWeights()
        DY2JetsFile = ROOT.TFile.Open(filelist["DY2JetsToLL"],"read")
        jetWeightMultiplicity["DY2JetsToLL"]=DY2JetsFile.Get("hWeights").GetSumOfWeights()
        DY3JetsFile = ROOT.TFile.Open(filelist["DY3JetsToLL"],"read")
        jetWeightMultiplicity["DY3JetsToLL"]=DY3JetsFile.Get("hWeights").GetSumOfWeights()
        DY4JetsFile = ROOT.TFile.Open(filelist["DY4JetsToLL"],"read")
        jetWeightMultiplicity["DY4JetsToLL"]=DY4JetsFile.Get("hWeights").GetSumOfWeights()

        W1JetsFile = ROOT.TFile.Open(filelist["W1JetsToLNu"],"read")
        jetWeightMultiplicity["W1JetsToLNu"]=W1JetsFile.Get("hWeights").GetSumOfWeights()

        W2JetsFile = ROOT.TFile.Open(filelist["W2JetsToLNu"],"read")
        jetWeightMultiplicity["W2JetsToLNu"]=W2JetsFile.Get("hWeights").GetSumOfWeights()
        W2JetsFileext1 = ROOT.TFile.Open(filelist["W2JetsToLNu_ext1"],"read")
        jetWeightMultiplicity["W2JetsToLNu"]+=W2JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W2JetsToLNu_ext1"]=W2JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W2JetsToLNu_ext1"]+=W2JetsFile.Get("hWeights").GetSumOfWeights()

        W3JetsFile = ROOT.TFile.Open(filelist["W3JetsToLNu"],"read")
        jetWeightMultiplicity["W3JetsToLNu"]=W3JetsFile.Get("hWeights").GetSumOfWeights()
        W3JetsFileext1 = ROOT.TFile.Open(filelist["W3JetsToLNu_ext1"],"read")
        jetWeightMultiplicity["W3JetsToLNu"]+=W3JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W3JetsToLNu_ext1"]=W3JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W3JetsToLNu_ext1"]+=W3JetsFile.Get("hWeights").GetSumOfWeights()

        W4JetsFile = ROOT.TFile.Open(filelist["W4JetsToLNu"],"read")
        W4JetsFileext1 = ROOT.TFile.Open(filelist["W4JetsToLNu_ext1"],"read")
        W4JetsFileext2 = ROOT.TFile.Open(filelist["W4JetsToLNu_ext2"],"read")
        jetWeightMultiplicity["W4JetsToLNu"]=W4JetsFile.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu"]+=W4JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu"]+=W4JetsFileext2.Get("hWeights").GetSumOfWeights()

        jetWeightMultiplicity["W4JetsToLNu_ext1"]=W4JetsFileext1.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu_ext1"]+=W4JetsFile.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu_ext1"]+=W4JetsFileext2.Get("hWeights").GetSumOfWeights()

        jetWeightMultiplicity["W4JetsToLNu_ext2"]=W4JetsFileext2.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu_ext2"]+=W4JetsFile.Get("hWeights").GetSumOfWeights()
        jetWeightMultiplicity["W4JetsToLNu_ext2"]+=W4JetsFileext1.Get("hWeights").GetSumOfWeights()


    Bkg = ["DY","W","TT","ST","EWK"]
    irBkg = ["ZZ","ZHToTauTau","vbf","WHTT"]
    TrialphaBkg = ["ttZ","ttW","WWZ","WZZ","ZZZ","WWW_4F","HZJ"]
    rareBkg = ["Other","rare","WZ"]
    finalDistributions = {}
    finalDistributions["Bkg"]=Bkg
    finalDistributions["data_obs"]=["data_obs"]
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
    if args.datameasureZH:
        finalDistributions["prompt1"]=["prompt1"]
        finalDistributions["prompt2"]=["prompt2"]
    if args.datadrivenZH or args.datadrivenSM:
        Bkg = ["FF"]
        rareBkg = ["DY","W","TT","ST","EWK","Other","rare","WZ"]
        finalDistributions["Bkg"]=Bkg
        finalDistributions["rareBkg"]=rareBkg

    allSkims = {}
    finalSkims ={}

    datadrivenPackage={}
    datadrivenPackage["bool"]=False


    if args.datadrivenZH or args.makeFakeHistos:
        print "MUST HAVE CREATED FF DISTRIBUTIONS BEFORE!"
        histodict = {}
        #inputFFile = ROOT.TFile.Open("FFskim_"+str(args.ffin)+".root","read")
        with uproot.open("skimmed_"+str(args.ffin)+".root") as inputFFile:
            histodict = createFakeFactorHistos(allcats, inputFFile)
        #inputFFile.Close()
        datadrivenPackage={}
        datadrivenPackage["bool"]=args.datadrivenZH
        fakemeasurefile = ROOT.TFile.Open("FFhistos_"+str(args.ffin)+"/fakemeasure.root","RECREATE")
        fakemeasurefile.cd()

        ss_1_tight = histodict[args.channel+"_FF_SS_1_tight"]["Nominal_data_obs"]["pt_3_ff"]
        ss_1_loose = histodict[args.channel+"_FF_SS_1_loose"]["Nominal_data_obs"]["pt_3_ff"]
        ss_2_tight = histodict[args.channel+"_FF_SS_2_tight"]["Nominal_data_obs"]["pt_4_ff"]
        ss_2_loose = histodict[args.channel+"_FF_SS_2_loose"]["Nominal_data_obs"]["pt_4_ff"]


        ss_1_tight_prompt = histodict[args.channel+"_FF_SS_1_tight"]["Nominal_prompt1"]["pt_3_ff"]
        ss_1_loose_prompt = histodict[args.channel+"_FF_SS_1_loose"]["Nominal_prompt1"]["pt_3_ff"]
        ss_2_tight_prompt = histodict[args.channel+"_FF_SS_2_tight"]["Nominal_prompt2"]["pt_4_ff"]
        ss_2_loose_prompt = histodict[args.channel+"_FF_SS_2_loose"]["Nominal_prompt2"]["pt_4_ff"]

        #subtracting prompt MC execpt for low stat channel
        # if args.channel!="mmem":
        ss_1_tight.Write("ss_1_tight",ROOT.TObject.kOverwrite)
        ss_1_loose.Write("ss_1_loose",ROOT.TObject.kOverwrite)
        ss_2_tight.Write("ss_2_tight",ROOT.TObject.kOverwrite)
        ss_2_loose.Write("ss_2_loose",ROOT.TObject.kOverwrite)


        ss_1_tight.Add(ss_1_tight_prompt,-1)
        ss_2_tight.Add(ss_2_tight_prompt,-1)
        ss_1_loose.Add(ss_1_loose_prompt,-1)
        ss_2_loose.Add(ss_2_loose_prompt,-1)




        f_1= ss_1_tight.Clone()
        f_1.Divide(ss_1_loose)
        f_2 = ss_2_tight.Clone()
        f_2.Divide(ss_2_loose)

        f_1.GetYaxis().SetTitleOffset(1.4)
        f_2.GetYaxis().SetTitleOffset(1.4)
        f_1.GetYaxis().SetMaxDigits(2)
        f_2.GetYaxis().SetMaxDigits(2)
        #ROOT.TGaxis().SetMaxDigits(2)

        f_1.SetName(args.channel+"_FakeRateLeg1")
        f_1.SetTitle(args.channel+" Fake Rate Measurement Leg1")
        f_1.GetXaxis().SetTitle("p_T Leg1")
        f_1.GetYaxis().SetTitle("Fake Rate for Leg1")
        f_2.SetName(args.channel+"_FakeRateLeg2")
        f_2.SetTitle(args.channel+" Fake Rate Measurement Leg2")
        f_2.GetXaxis().SetTitle("p_T Leg2")
        f_2.GetYaxis().SetTitle("Fake Rate for Leg2")

        tf_1 = ROOT.TF1("tf_1","[0]",f_1.GetXaxis().GetXmin(),f_1.GetXaxis().GetXmax())
        tf_2 = ROOT.TF1("tf_2","[0]",f_2.GetXaxis().GetXmin(),f_2.GetXaxis().GetXmax())

        c=ROOT.TCanvas("canvas","",0,0,600,600)
        #ROOT.gStyle.SetOptFit()
        ROOT.gStyle.SetOptStat(0)
        #ROOT.gStyle.SetOptFit(1)
        f_1.Draw()
        f_1.Fit("tf_1")
        c.SaveAs("FFhistos_"+str(args.ffin)+"/"+args.channel+"_fakerate1.png")
        f_2.Draw()
        f_2.Fit("tf_2")
        c.SaveAs("FFhistos_"+str(args.ffin)+"/"+args.channel+"_fakerate2.png")
        tf_1.Write(tf_1.GetName(),ROOT.TObject.kOverwrite)
        tf_2.Write(tf_2.GetName(),ROOT.TObject.kOverwrite)
        f_1.Write(f_1.GetName(),ROOT.TObject.kOverwrite)
        f_2.Write(f_2.GetName(),ROOT.TObject.kOverwrite)

        ss_1_tight.Write("ss_1_tight_sub",ROOT.TObject.kOverwrite)
        ss_1_loose.Write("ss_1_loose_sub",ROOT.TObject.kOverwrite)
        ss_2_tight.Write("ss_2_tight_sub",ROOT.TObject.kOverwrite)
        ss_2_loose.Write("ss_2_loose_sub",ROOT.TObject.kOverwrite)
        ss_1_tight_prompt.Write("ss_1_tight_prompt",ROOT.TObject.kOverwrite)
        ss_1_loose_prompt.Write("ss_1_loose_prompt",ROOT.TObject.kOverwrite)
        ss_2_tight_prompt.Write("ss_2_tight_prompt",ROOT.TObject.kOverwrite)
        ss_2_loose_prompt.Write("ss_2_loose_prompt",ROOT.TObject.kOverwrite)




        datadrivenPackage["fakerate1"]=f_1.Clone()
        datadrivenPackage["fakerate2"]=f_2.Clone()
        datadrivenPackage["fitrate1"]=tf_1
        datadrivenPackage["fitrate2"]=tf_2
        datadrivenPackage["fakemeasurefile"]=fakemeasurefile
        #fakemeasurefile.Close()

    #EventWeights = getEventWeightDicitonary()
    #importing special fake factor class from parametrization
    fakefactorObj= fakefactor()
    if args.datadrivenSM:
        #fakefactorObj = fakefactor("/eos/home-s/shigginb/fakefactor/")
        #fakefactorObj.loadHistograms("/eos/home-s/shigginb/fakefactors/",args.year)
        fakefactorObj.loadHistograms("./fakefactors/",args.year)
        #fakefactorObj.Print()

    #exit()
    # ROOT.fail
    return allcats, HAA_processes,finalDistributions,weightHistoDict,jetWeightMultiplicity,datadrivenPackage,fakefactorObj

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
       print "events tree not present"
   finalSkims = {s:{c: dict() for c in allcats.keys()} for s in systematics}


   rootFiles={sys:glob.glob(rootfiledir+"/"+sys+"_*") for sys in systematics}


   mainTList = {}
   mainOutputTree={}
   emptyFiles=[]
   nonemptyFiles=[]
   for sys, globfiles in rootFiles.iteritems():
      for globfile in globfiles:
         print "working on glob ",globfile,"  size   ",os.path.getsize(globfile) 
         for cat,catObj in allcats.iteritems():
            for nickname, processObj in processes.iteritems():
                for process in processObj.cuts.keys():
                    if str(globfile).split("/")[1]==sys+"_"+cat+"_"+nickname+"_"+process:
                        print "found match ",sys+"_"+cat+"_"+nickname+"_"+process
                        fin = uproot.open(globfile)
                        if len(fin[sys+"_"+cat+"_"+nickname+"_"+process].keys())==12: nonemptyFiles.append(sys+"_"+cat+"_"+nickname+"_"+process)
                        else: emptyFiles.append(sys+"_"+cat+"_"+nickname+"_"+process)
                            
                        

                        #with uproot.open(globfile) as fin:
                        #    tree = fin[sys+"_"+cat+"_"+nickname+"_"+process]
                        #    mainArrays = tree.arrays()
                        #    # mainArrays={}
                        #    # for branch,data in tree.iteritems():
                        #    #     mainArrays[branch]=data.array(library="np")
                        #    #mainArrays["nickname"]=mainArrays["nickname"].array(library="np")
                        #    #print "tree ",sys+"_"+process," entries ",len(mainArrays["mll"])
                        #    for catDist, final in finalDistributions.iteritems():
                        #        for processOut in final:
                        #            if (processOut==process) and (catDist not in finalSkims[sys][cat]):
                        #                #print "first output for process ",process," finalDist cat ",catDist
                        #                finalSkims[sys][cat][catDist] = mainArrays
                        #                continue
                        #            elif (processOut==process) and (catDist in finalSkims[sys][cat]):
                        #                #print "adding to finalskims ", catDist,"  for process ",process," finalDist cat ",catDist
                        #                for branch in finalSkims[sys][cat][catDist].keys():
                        #                    if branch=="nickname":
                        #                         print "nickname branch first ele ",finalSkims[sys][cat][catDist][branch][0]
                        #                         print "nickname branch whole     ",finalSkims[sys][cat][catDist][branch]
                        #                        #finalSkims[sys][cat][catDist][branch]=finalSkims[sys][cat][catDist][branch].astype("S40")
                        #                    try:
                        #                        #finalSkims[sys][cat][catDist][branch]=np.concatenate((finalSkims[sys][cat][catDist][branch],mainArrays[branch]))
                        #                    except:
                        #                        print "possible missing branch ",branch
                        #            else:
                        #                continue

   #skimFile = ROOT.TFile("skimmed_"+outputstring+".root","recreate")
   #skimFile.cd()
   #for cat, catObj in allcats.iteritems():
   #    skimFile.cd()
   #    skimFile.mkdir(cat)
   #    skimFile.cd(cat)
   #    for sys in finalSkims.keys():
   #        dataTypes =[[],[]]
   #        random_sample = finalSkims[sys][cat].values()[0]
   #        for branch in random_sample.keys():
   #            if branch=="nickname":
   #                branchdatatype = "S40"
   #            else:
   #                branchdatatype = random_sample[branch].dtype
   #            dataTypes[0].append(branch)
   #            dataTypes[1].append(branchdatatype)
   #        for catDist in finalSkims[sys][cat].keys():
   #            data = np.zeros(len(finalSkims[sys][cat][catDist][branch]),dtype={'names':dataTypes[0],'formats':dataTypes[1]})
   #            for branch in data.dtype.names:
   #                if len(finalSkims[sys][cat][catDist][branch].shape) == 1:   # flat array important for jagged arrays of input data
   #                    data[branch] = finalSkims[sys][cat][catDist][branch]
   #                else:
   #                    data[branch] = finalSkims[sys][cat][catDist][branch][:,0]
   #            treeOut = root_numpy.array2tree(data, name=sys+"_"+catDist)
   #            treeOut.Write()
   #skimFile.Close()

   return nonemptyFiles,emptyFiles

if __name__ == "__main__":
    import datetime
    begin_time = datetime.datetime.now()
    from multiprocessing import *
    import multiprocessing as mp
    import logging
    import os
    import shutil


    from utils.Weights import jet_exc_samples


    import argparse


    parser = argparse.ArgumentParser(description="This file generates root files containing Histograms ... files in utils contain selections and settings")
    parser.add_argument("-o",  "--outname", default="",  help="postfix string")
    parser.add_argument("-fi",  "--ffin", default="",  help="fake factor files")
    parser.add_argument("-year",  "--year", default="2016",  help="Year")
    parser.add_argument("-fo",  "--ffout", default="",  help="fake factor files to output")
    parser.add_argument("-c",  "--categories", default="categories_array.yaml",  help="categories yaml file")
    parser.add_argument("-ch",  "--channel", default="mmmt",  help="Please list the channel for fake factor histograms")
    parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
    parser.add_argument("-i",  "--dir", default="/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov7_basic_10_6_4/src/2016_v7/",  help="Input files")
    parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
    parser.add_argument("-dm",  "--datameasure", default=False,action='store_true',  help="Use DataDriven Method measure part")
    parser.add_argument("-dmZH",  "--datameasureZH", default=False,action='store_true',  help="Use DataDriven Method measure part")
    parser.add_argument("-dbg",  "--debug", default=False,action='store_true',  help="Disable the parallel processing mode")
    parser.add_argument("-sys",  "--systematics", default=False,action='store_true',  help="Run on systematic trees, don't run nominal at the same time")
    parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-ddSM",  "--datadrivenSM", default=False,action='store_true',  help="Use DataDriven Method")
    parser.add_argument("-ex",  "--extract", default=False,action='store_true',  help="Additional Cuts for Extraction")
    parser.add_argument("-ff",  "--makeFakeHistos", default=False,action='store_true',  help="Just make fake rate histos")
    parser.add_argument("-v",  "--verbose", default=False,action='store_true',  help="print per event")
    parser.add_argument("-t",  "--test", default=False,action='store_true',  help="only do 1 event to test code")
    parser.add_argument("-mt",  "--mt", default=False,action='store_true',  help="Use Multithreading")
    parser.add_argument("-pt",  "--maxprint", default=False,action='store_true',  help="Print Info on cats and processes")
    parser.add_argument("-combine",  "--combine", default=False,action='store_true',  help="just combine root files")
    parser.add_argument("-pc",  "--processingcores", default=12,type=int,help="Number of cores for multiprocessing")

    args = parser.parse_args()
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
        systematics =[ "scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
                       "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
                       "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
                       "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"]
    else:
        systematics =[ "Events"]
    #systematics =[ "Events"]

    jet_exc_samples = jet_exc_samples[args.year]

    print "combining root files"

    nonemptyFiles,emptyFiles=\
    checkRootFiles(systematics, allcats, HAA_processes,\
                     finalDistributions, "massOutputDir_"+args.outname,\
                     args.channel, args.ffout)
    print "emptyFiles ",emptyFiles
    #print "nonemptyFiles ",nonemptyFiles

    try:
        datadrivenPackage["fakemeasurefile"].Close()
    except:
        print "problem closing fake measurement helper file ... shutting down anyway"

    print("computation time")
    print(datetime.datetime.now() - begin_time)
    print("arguments used")
    print(args)
