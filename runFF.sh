echo "Running Fake Factor Sequence ..."
echo "Applying Fake Factors ..."
python MakeDistributions_HAA_2016.py -c cat_mmmt_2016.yaml -csv MCsamples_2016_v6_yaml.csv -i /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/ -p processes_special_mmmt.yaml -o 2016_skimff_mmmt -ch mmmt -s -ddZH -fi 2016_ff_mmmt

echo "Making Plots ..."
python MakePlots_skimmed.py -i skimmed_2016_skimff_mmmt.root -o 2016_skimmed_ff -c cat_mmmt_2016.yaml --ch mmmt -ddZH



echo "Copying Plots Over ..."
cp -r outplots_2016_skimmed_ff/ /eos/home-s/shigginb/HAA_Plots/.
