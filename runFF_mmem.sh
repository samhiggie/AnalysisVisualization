echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
python MakeDataCards_array_HAA.py -c cat_mmem_2016.yaml -csv MCsamples_2016_v6_yaml.csv -i /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/ -p processes_special_mmem.yaml -dmZH -o 2016_dm_mmem -fo 2016_mmem -ch mmem

echo "Copying over Fake Factors ..."
cp -r out2016_dm_mmem/pt_*.root FFhistos_2016_mmem/.

echo "Applying Fake Factors ..."
python MakeDataCards_array_HAA.py -c cat_mmem_2016.yaml -csv MCsamples_2016_v6_yaml.csv -i /afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/ -p processes_special_mmem.yaml -ddZH -o 2016_mmem -fi 2016_mmem -ch mmem

echo "Making Plots ..."
python MakePlots_histos.py -i 2016_mmem -o 2016_mmem -ddZH -mhs -c cat_mmem_2016.yaml --ch mmem


echo "Copying Plots Over ..."
cp FFhistos_2016_mmem/*.png outplots_2016_mmem_mmem_inclusive/.

echo "Copying fake rate plots Over ..."
cp -r outplots_2016_mmem_mmem_inclusive /eos/home-s/shigginb/HAA_Plots/.
