echo "Running Fake Factor Sequence ..."
echo "Measuring Fake Factors ..."
python MakeDataCards_HAA.py -o n6_datadriven_measure_ZH_loose -dmZH

echo "Copying over Fake Factors ..."
cp outn6_datadriven_measure_ZH_loose/*pt_*.root FFhistos/.

echo "Applying Fake Factors ..."
python MakeDataCards_HAA.py -o n6_datadriven_ZH_loose -ddZH

echo "Making Plots ..."
python MakePlots.py -i n6_datadriven_ZH_loose -o n6_datadriven_ZH_loose -ddZH


echo "Copying Plots Over ..."
cp -rf outplots_n6_datadriven_ZH_loose_mmmt_inclusive/ /eos/home-s/shigginb/HAA_Plots/.

echo "Copying fake rate plots Over ..."
cp FFhistos/fakerate*.png /eos/home-s/shigginb/HAA_Plots/outplots_n6_datadriven_ZH_loose_mmmt_inclusive/.
