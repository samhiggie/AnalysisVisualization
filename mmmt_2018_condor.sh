#!/bin/bash
cd ${_CONDOR_SCRATCH_DIR}
export X509_USER_PROXY=$1
voms-proxy-info -all
voms-proxy-info -all -file $1
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
eval `scramv1 project CMSSW CMSSW_10_2_16_patch1`
cd CMSSW_10_2_16_patch1/src
eval `scramv1 runtime -sh`
git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs
cd ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/
cp ${_CONDOR_SCRATCH_DIR}/* .
mkdir ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/utils
mkdir ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/
cp __init__.py ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/utils/.
cp Parametrization.py ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/utils/.
cp Weights.py ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/utils/.
cp functions.py ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/utils/.
cp FitHistograms_eleFR_2016.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_eleFR_2016_sub09.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_eleFR_2016_sub11.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_eleFR_2017.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_eleFR_2017_sub09.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_eleFR_2017_sub11.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_eleFR_2018.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_eleFR_2018_sub09.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_eleFR_2018_sub11.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_FR.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_muFR_2016.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_muFR_2016_sub09.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_muFR_2016_sub11.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_muFR_2017.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_muFR_2017_sub09.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_muFR_2017_sub11.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_muFR_2018.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_muFR_2018_sub09.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_muFR_2018_sub11.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_tauFR_2016.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_tauFR_2017.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
cp FitHistograms_tauFR_2018.root ${_CONDOR_SCRATCH_DIR}/CMSSW_10_2_16_patch1/src/fakefactors/.
scram b -j 4
eval `scramv1 runtime -sh`
ls -altrh
channel=mmmt
year=2018
input=/eos/home-s/shigginb/HAA_ntuples/2018_paper/
process=processes_special_mmmt.yaml
output1=2018_condortry_mmmt
output1ex=2018_condortryex_mmmt
fo=2018_condortry_dm_mmmt
csv=MCsamples_2018_paper.csv
cat=cat_mmmt_2018.yaml
catsys=cat_sys_mmmt_2018.yaml
python MakeDistributions_v7.py -c $catsys -csv $csv -i $input -p $process -o ${output1}_sys -ch ${channel} -s -sys -fi $fo -fo ${output1}_sys -year ${year} -ex
python MakeDistributions_v7.py -c $catsys -csv $csv -i $input -p $process -o $output1 -ch ${channel} -s -ddSM -fi $fo -fo $output1 -year $year
hadd -f skimmed_${output1}_combined.root skimmed_${output1}.root skimmed_${output1}_sys.root
python MakePlots_skimmed_sys.py -i skimmed_${output1}_combined.root -o $output1 -c $catsys --ch ${channel} -p $process -de -ddZH -year $year
python MakeDistributions_v7.py -c $catsys -csv $csv -i $input -p $process -o ${output1ex}_sys -ch ${channel} -s -sys -fi $fo -fo ${output1ex}_sys -year ${year} -ex
python MakeDistributions_v7.py -c $catsys -csv $csv -i $input -p $process -o $output1ex -ch ${channel} -s -ddSM -fi $fo -fo $output1ex -year $year -ex
hadd -f skimmed_${output1ex}_combined.root skimmed_${output1ex}.root skimmed_${output1ex}_sys.root
python MakePlots_skimmed_sys.py -i skimmed_${output1ex}_combined.root -o $output1ex -c $catsys --ch ${channel} -p $process -de -ddZH -year $year
cp -r outplots_${output1}_*/ /eos/home-s/shigginb/HAA_Plots/.
cp -r outplots_${output1ex}_*/ /eos/home-s/shigginb/HAA_Plots/.
cp  skimmed_${output1ex}_combined.root ${_CONDOR_SCRATCH_DIR}/.
cp -r outplots_${output1}_*/   ${_CONDOR_SCRATCH_DIR}/.
cp -r outplots_${output1ex}_*/ ${_CONDOR_SCRATCH_DIR}/.
