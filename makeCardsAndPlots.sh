echo "Doing it alllll ! "
python MakeDataCards_HAA.py
python MakePlots.py
#outstuff=/eos/home-s/shigginb/HAA_Plots/Inc/.
#echo "copy-ing files over to" $outstuff
#cp outplots/*.png $outstuff
cp -r outplots_mmmt_* /eos/home-s/shigginb/HAA_Plots/.
