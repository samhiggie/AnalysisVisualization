echo "Running Fake Factor Sequence ..."
echo "Applying Fake Factors ..."
python MakeDistributions_HAA_2016.py -c cat_mmmt_2016.yaml -csv MCsamples_2016_v6_yaml.csv -i /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/ -p processes_special_mmmt.yaml -o 2016_skim_mc_mmmt -ch mmmt -s -fi 2016_ff_mmmt

echo "Making Plots ..."
python MakePlots_skimmed.py -i skimmed_2016_skim_mc_mmmt.root -o 2016_skimmed_mc -c cat_mmmt_2016.yaml --ch mmmt



echo "Copying Plots Over ..."
cp -r outplots_2016_skimmed_mc/ /eos/home-s/shigginb/HAA_Plots/.
