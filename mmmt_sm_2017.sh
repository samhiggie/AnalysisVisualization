echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
channel=mmmt
year=2017

#input=/eos/home-s/shigginb/HAA_ntuples/${year}_v7/
input=/eos/home-s/shigginb/HAA_ntuples/${year}_paper/


process=processes_trigger_${channel}.yaml
#csv=MCsamples_${year}_v7.csv
#csv=MCsamples_${year}_v7_final.csv
#csv=samples_${year}_v7.csv
csv=MCsamples_${year}_paper.csv

echo "using csv "$csv
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
fo=${year}_${mainout}_dm_${channel}
cat=cat_${channel}_${year}.yaml

catsys=cat_sys_${channel}_${year}.yaml


#python MakeDistributions_v7.py -c $catsys -csv $csv -i $input -p $process -o ${output1}_sys -ch ${channel} -s -sys -fi $fo -fo ${output1}_sys -year ${year}


#exit 1
#measure
#python MakeDistributions_v7.py -c $cat  -csv $csv  -i $input -p $process -dmZH -o $output0 -fo $fo -ch ${channel} -year $year
#python MakeDistributions_v7.py -c $cat  -csv $csv  -i $input -p $process -o $output0 -fo $fo -ch ${channel}
#python MakePlots_skimmed_sys.py -i skimmed_${output0}.root -o $output0 -c $cat --ch ${channel} -p $process

echo "Copying over Fake Factors ..."
cp -r out$output0/pt_*.root FFhistos_$fo/.

#echo "Applying Fake Factors ..."
#python MakeDistributions_v7.py -c $cat -csv $csv -i $input -p $process -o ${output1} -ch ${channel} -s -ddZH -fi $fo -fo $output1  -year $year
python MakeDistributions_v7.py -c $cat -csv $csv -i $input -p $process -o ${output1} -ch ${channel} -s -ddSM -fi $fo -fo $output1  -year $year
#extraction
#python MakeDistributions_v7.py -c $cat -csv $csv -i $input -p $process -o ${output1} -ch ${channel} -s -ddZH -fi $fo -fo $output1  -year $year -ex

#echo "Making Plots ..."
python MakePlots_skimmed_sys.py -i skimmed_${output1}.root -o ${output1} -c $cat --ch ${channel} -p $process -de -ddZH -year $year


echo "combinging systematics and nominal distributions"
#hadd skimmed_${output1}_combined.root skimmed_${output1}.root skimmed_${output1}_sys.root

#echo "Copying Plots Over ..."
cp -r FFhistos_$fo/*.png outplots_${output1}_Nominal/.
cp -r outplots_${output1}_Nominal/ /eos/home-s/shigginb/HAA_Plots/.


#cp -r outplots_${output0}_Nominal/ /eos/home-s/shigginb/HAA_Plots/.
#cp -r outplots_${output1}_cats/ /eos/home-s/shigginb/HAA_Plots/.

rm -rf massOutputDir_${year}_${mainout}_dm_${channel}
rm -rf massOutputDir_${year}_${mainout}_${channel}_sys
rm -rf massOutputDir_${year}_${mainout}_${channel}
