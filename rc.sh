#/bin/bash.sh
#for channel in {mmtt,mmet,mmem,mmmt}; do echo ${channel}; hadd -f skimmed_comboex_${channel}.root skimmed_201*_comboex_${channel}_combined.root; done;

#for channel in {mmtt,mmet,mmem,mmmt}; do echo ${channel}; hadd -f skimmed_combo_${channel}.root skimmed_201*_combo_${channel}_combined.root; done;

for channel in {mmtt,mmet,mmem,mmmt};
do echo ${channel};
python MakePlots_skimmed_sys.py -i skimmed_combo_${channel}.root -o ${channel}_combo -c cat_sys_${channel}_2018.yaml --ch ${channel} -p processes_special_${channel}.yaml -de -ddZH -year RunII -co;
done;

#hadd -f skimmed_combo.root skimmed_combo_mmtt.root skimmed_combo_mmet.root skimmed_combo_mmmt.root skimmed_combo_mmem.root

#hadd -f skimmed_comboex.root skimmed_comboex_mmtt.root skimmed_comboex_mmet.root skimmed_comboex_mmmt.root skimmed_comboex_mmem.root

python MakePlots_skimmed_sys.py -i skimmed_combo.root -o allyears_combo -c cat_sys_mmmt_2018.yaml --ch mmmt -p processes_special_mmmt.yaml -de -ddZH -year RunII -co -ac

cp -r outplots_mm*_combo_Nominal/ /eos/home-s/shigginb/HAA_Plots/.
cp -r outplots_allyears_combo_Nominal/ /eos/home-s/shigginb/HAA_Plots/.
