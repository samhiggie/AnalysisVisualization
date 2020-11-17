echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
#common
input=/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/
process=processes_special_mmmt.yaml
csv=MCsamples_2016_v6_yaml.csv
#channel specific
output0=2016_prompt_dm_mmmt
output1=2016_prompt_mmmt
fo=2016_prompt_mmmt
cat=cat_mmmt_2016.yaml

python MakeDistributions_HAA_2016.py -c $cat  -csv $csv  -i $input -p $process -dmZH -o $output0 -fo $fo -ch mmmt

echo "Copying over Fake Factors ..."
cp -r out$output0/pt_*.root FFhistos_$fo/.

echo "Applying Fake Factors ..."
python MakeDistributions_HAA_2016.py -c $cat -csv $csv -i $input -p $process -o $output1 -ch mmmt -s -ddZH -fi $fo

echo "Making Plots ..."
python MakePlots_skimmed.py -i skimmed_$output1.root -o $output1 -c $cat --ch mmmt -ddZH

echo "Copying Plots Over ..."
cp -r FFhistos_$fo/*.png outplots_$output1/.
cp -r outplots_$output1/ /eos/home-s/shigginb/HAA_Plots/.
