echo "Installing required dependencies ..."

cd ${CMSSW_BASE}/src

git clone https://github.com/CMS-HTT/RecoilCorrections.git HTT-utilities/RecoilCorrections

git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs

scram b -j 21

cmsenv

cd ${CMSSW_BASE/src/AnalysisVisualization}
