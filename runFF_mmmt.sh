echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
python MakeDataCards_array_HAA.py -c cat_mmmt_2016.yaml -csv MCsamples_2016_v6_yaml.csv -i /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/ -p processes_special.yaml -dmZH -o 2016_ff_dm_mmmt -fo 2016_ff_mmmt -ch mmmt

echo "Copying over Fake Factors ..."
cp -r out2016_ff_dm_mmmt/pt_*.root FFhistos_2016_ff_mmmt/.

echo "Applying Fake Factors ..."
python MakeDataCards_array_HAA.py -c cat_mmmt_2016.yaml -csv MCsamples_2016_v6_yaml.csv -i /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/ -p processes_special.yaml -ddZH -o 2016_ff_mmmt -fi 2016_ff_mmmt -ch mmmt

echo "Making Plots ..."
python MakePlots_histos.py -i 2016_ff_mmmt -o 2016_ff_mmmt -ddZH -mhs -c cat_mmmt_2016.yaml --ch mmmt


echo "Copying Plots Over ..."
cp FFhistos_2016_ff_mmmt/*.png outplots_2016_ff_mmmt_mmmt_inclusive/.

echo "Copying fake rate plots Over ..."
cp -r outplots_2016_ff_mmmt_mmmt_inclusive /eos/home-s/shigginb/HAA_Plots/.
