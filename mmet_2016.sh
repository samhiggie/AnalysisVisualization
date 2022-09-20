channel=mmet
year=2016

input=/eos/home-s/shigginb/HAA_ntuples/${year}_paper/


process=processes_special_${channel}.yaml
csv=MCsamples_${year}_paper.csv
if [ -z $1 ]
then
    mainout='test'
    echo "User may want to provide an output string via single argument"
else
    mainout=$1
    echo 'using '$mainout' as output'
fi
output0=${year}_${mainout}_dm_${channel}
output1=${year}_${mainout}_${channel}
output1ex=${year}_${mainout}ex_${channel}
fo=${year}_${mainout}_dm_${channel}
fmc=${year}_${mainout}_${channel}
cat=cat_${channel}_${year}.yaml

cat=cat_sys_${channel}_${year}.yaml
catsys=cat_sys_${channel}_${year}.yaml


python MakeDistributions_v7.py -c $catsys -csv $csv -i $input -p $process -o ${output1}_sys -ch ${channel} -s -sys -fi $fo -fo ${output1}_sys -year ${year}

python MakeDistributions_v7.py -c $catsys -csv $csv -i $input -p $process -o $output1 -ch ${channel} -s -ddSM -fi $fo -fo $output1 -year $year

hadd -f skimmed_${output1}_combined.root skimmed_${output1}.root skimmed_${output1}_sys.root

python MakePlots_skimmed_sys.py -i skimmed_${output1}_combined.root -o $output1 -c $catsys --ch ${channel} -p $process -de -ddZH -year $year




python MakeDistributions_v7.py -c $catsys -csv $csv -i $input -p $process -o ${output1ex}_sys -ch ${channel} -s -sys -fi $fo -fo ${output1ex}_sys -year ${year} -ex


python MakeDistributions_v7.py -c $catsys -csv $csv -i $input -p $process -o $output1ex -ch ${channel} -s -ddSM -fi $fo -fo $output1ex -year $year -ex

hadd -f skimmed_${output1ex}_combined.root skimmed_${output1ex}.root skimmed_${output1ex}_sys.root

python MakePlots_skimmed_sys.py -i skimmed_${output1ex}_combined.root -o $output1ex -c $catsys --ch ${channel} -p $process -de -ddZH -year $year




echo "Copying Plots Over ..."
cp -r outplots_${output1}_*/ /eos/home-s/shigginb/HAA_Plots/.
cp -r outplots_${output1ex}_*/ /eos/home-s/shigginb/HAA_Plots/.

rm -rf massOutputDir_${year}_${mainout}_dm_${channel}
rm -rf massOutputDir_${year}_${mainout}_${channel}_sys
rm -rf massOutputDir_${year}_${mainout}ex_${channel}
rm -rf massOutputDir_${year}_${mainout}ex_${channel}_sys
rm -rf massOutputDir_${year}_${mainout}_${channel}
