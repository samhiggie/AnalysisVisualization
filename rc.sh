#/bin/bash.sh
#for channel in {mmtt,mmet,mmem}; do echo ${channel}; hadd -f skimmed_combo_${channel}.root skimmed_201*_combo_${channel}_combined.root; done;

for channel in {mmtt,mmet,mmem};
do echo ${channel};
python MakePlots_skimmed_sys.py -i skimmed_combo_${channel}.root -o ${channel}_combo -c cat_sys_${channel}_2018.yaml --ch ${channel} -p processes_special_${channel}.yaml -de -ddZH -year RunII -co;
done;

hadd -f skimmed_combo.root skimmed_combo_*.root

python MakePlots_skimmed_sys.py -i skimmed_combo.root -o allyears_combo -c cat_sys_mmmt_2018.yaml --ch mmmt -p processes_special_mmmt.yaml -de -ddZH -year RunII -co
