echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
#common
input=/afs/cern.ch/work/s/shigginb/cmssw/HAA/nanov6_10_2_9/src/nano6_2016/
process=processes_special_mmmt.yaml
csv=MCsamples_2016_v6_yaml.csv
#channel specific
output0=2016_prompt_dm_mmmt
output1=2016_full_mmmt
output2=2016_full_mmmt_histos
fo=2016_prompt_mmmt
cat=cat_mmmt_2016.yaml
inc=_mmmt_inclusive
f1t=_mmmt_FF_SS_1_tight
f1l=_mmmt_FF_SS_1_loose
f2t=_mmmt_FF_SS_2_tight
f2l=_mmmt_FF_SS_2_loose

#python MakeDistributions_HAA_2016.py -c $cat  -csv $csv  -i $input -p $process -dmZH -o $output0 -fo $fo -ch mmmt
#
#echo "Copying over Fake Factors ..."
#cp -r out$output0/pt_*.root FFhistos_$fo/.
#
echo "Applying Fake Factors ..."
python MakeDistributions_HAA_2016.py -c $cat -csv $csv -i $input -p $process -o $output1 -ch mmmt -s -ddZH -fi $fo

#python MakeDistributions_HAA_2016.py -c $cat  -csv $csv  -i $input -p $process -ddZH -o $output1 -fi $fo -ch mmmt

#echo "Making Plots ..."
#python MakePlots_skimmed.py -i skimmed_$output1$inc.root -o $output1$inc -c $cat --ch mmmt -ddZH

echo "Making FF control plots"
#python MakePlots_histos.py -i $output1 -o $output2 -c $cat --ch mmmt -ddZH
python MakePlots_skimmed_cats.py -i $output1 -o $output1 -c $cat --ch mmmt -ddZH

echo "Copying Plots Over ..."
cp -r FFhistos_$fo/*.png outplots_$output1/.
cp -r outplots_$output1/ /eos/home-s/shigginb/HAA_Plots/.
