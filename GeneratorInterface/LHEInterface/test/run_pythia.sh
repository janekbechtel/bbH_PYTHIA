for SEED in 01 02 03 04 05 06 07 08 09 10 11 12
do
cmsRun POWHEG_Pythia8_bbh_tautau_cfg.py $SEED & 
done
wait

hadd -f -n 12 bbh_mumu_700gev_powheg.root bbh_mumu_700gev_powheg_*.root


