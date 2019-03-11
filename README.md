# bbH_studies

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch

source $VO_CMS_SW_DIR/cmsset_default.sh

scram project CMSSW_7_1_23

cd CMSSW_7_1_23/src

cmsenv

git clone https://github.com/janekbechtel/bbH_studies .

scramv1 b -j 12

cd GeneratorInterface/LHEInterface/test/

cmsRun POWHEG_Pythia8_bbh_tautau_cfg.py

#OR
#cmsRun aMCNLO_MG5_Pythia8_bbh_tautau_cfg.py

