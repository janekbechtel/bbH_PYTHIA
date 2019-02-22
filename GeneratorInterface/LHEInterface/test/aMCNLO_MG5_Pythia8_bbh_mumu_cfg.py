#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")

process.source = cms.Source("LHESource",
# mH = 300 GeV, PDF4LHC15_nf4_30
#	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_3_0_beta/bbH_4FS_yb2/Events/run_04/events.lhe')
# mH = 500 GeV, PDF4LHC15_nf4_30
#	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_3_0_beta/bbH_4FS_yb2/Events/run_05/events.lhe')
# mH = 700 GeV, PDF4LHC15_nf4_30
#	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_3_0_beta/bbH_4FS_yb2/Events/run_06/events.lhe')
# mH = 700 GeV, PDF4LHC15_nf4_30, alpha x sqrt(2)
#	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_3_0_beta/bbH_4FS_yb2/Events/run_11/events.lhe')
# mH = 700 GeV, PDF4LHC15_nf4_30, alpha / sqrt(2)
#	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_3_0_beta/bbH_4FS_yb2/Events/run_13/events.lhe')
#	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_6_1/CPggH/Events/run_01/events.lhe')
#	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_6_1/CPggH/Events/run_17/events.lhe')
#
#
#
# Marius events with slpha 0.25
	fileNames = cms.untracked.vstring('file:/localscratch/bba_mumu_lhe/bba_mumu_30gev_yb2.lhe')
#	fileNames = cms.untracked.vstring('file:/localscratch/bba_mumu_lhe/bba_mumu_30gev_ybyt.lhe')
#
# POWHEG SM gg->h with t+b
#  fileNames = cms.untracked.vstring('file:/localscratch/POWHEG-BOX/gg_H_quark-mass-effects/testrun-lhc/sm_gg_h_mumu_30gev_powheg.lhe')
#    fileNames = cms.untracked.vstring('file:/localscratch/POWHEG-BOX/gg_H_MSSM/testrun-lhc-A/mssm_gg_a_mumu_ma30tanb40_powheg_hfact30gev.lhe')
#     fileNames = cms.untracked.vstring('file:/localscratch/POWHEG-BOX/gg_H_MSSM/testrun-lhc-A/mssm_gg_a_mumu_ma30tanb40_powheg_hfact15gev.lhe')
#     fileNames = cms.untracked.vstring('file:/localscratch/POWHEG-BOX/gg_H_MSSM/testrun-lhc-A/mssm_gg_a_mumu_ma30tanb40_powheg_hfact7gev.lhe')
#      fileNames = cms.untracked.vstring('file:/localscratch/POWHEG-BOX/gg_H_MSSM/testrun-lhc-A/mg5-amcnlo-SM_ggA_with_b_loop_only_marius_hfact15gev.lhe')
#      fileNames = cms.untracked.vstring('file:/localscratch/POWHEG-BOX/gg_H_MSSM/testrun-lhc-A/mg5-amcnlo-SM_ggA_b_loop_g_bb_oliver-mattelaer_v2.lhe')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.configurationMetadata = cms.untracked.PSet(
	version = cms.untracked.string('alpha'),
	name = cms.untracked.string('LHEF input'),
	annotation = cms.untracked.string('ttbar')
)

process.load("Configuration.StandardSequences.Services_cff")

process.RandomNumberGeneratorService.generator = cms.PSet(
	initialSeed = cms.untracked.uint32(123456789),
	engineName = cms.untracked.string('HepJamesRandom')
)


from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

# gen jets ak4
#
process.load("RecoJets.JetProducers.ak4GenJets_cfi")
process.ak4Jets  =  process.ak4GenJets.clone()
process.ak4Jets.src = cms.InputTag("genParticlesForAK4Jets")
#
#input particles for gen jets 
#
process.load("RecoJets.Configuration.GenJetParticles_cff")

process.genParticlesForAK4Jets = process.genParticlesForJets.clone()
process.genParticlesForAK4Jets.excludeFromResonancePids = cms.vuint32(11,12,13,14,15,16)
process.genParticlesForAK4Jets.ignoreParticleIDs +=       cms.vuint32(11,12,13,14,15,16)

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
	eventsToPrint = cms.untracked.uint32(2),
        UseExternalGenerators = cms.untracked.bool(True),
	maxEventsToPrint = cms.untracked.int32(2),
    	pythiaPylistVerbosity = cms.untracked.int32(1),
    	filterEfficiency = cms.untracked.double(1.0),
    	pythiaHepMCVerbosity = cms.untracked.bool(False), 
        comEnergy = cms.double(13000.),
	SLHAFileForPythia8 = cms.string('powheg-fh-output_mumu.slha'),
    PythiaParameters = cms.PSet(
        pythia8CommonSettings = cms.vstring(
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = on',
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 10.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on'),
        pythia8CUEP8M1Settings = cms.vstring(
            'Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:pT0Ref=2.4024', 
            'MultipartonInteractions:ecmPow=0.25208', 
            'MultipartonInteractions:expPow=1.6'),
        processParameters = cms.vstring(
	    'SpaceShower:pTmaxMatch = 1',        # aMC@NLO-PY8 interface settings
            'SpaceShower:pTmaxFudge = 1',        # aMC@NLO-PY8 interface settings
	    'TimeShower:pTmaxMatch = 1',         # aMC@NLO-PY8 interface settings
	    'TimeShower:pTmaxFudge = 1',         # aMC@NLO-PY8 interface settings
	    'SpaceShower:MEcorrections = off',   # aMC@NLO-PY8 interface settings
	    'TimeShower:MEcorrections = off',    # aMC@NLO-PY8 interface settings
	    'TimeShower:globalRecoil = on' ,     # aMC@NLO-PY8 interface settings
            'TimeShower:weightGluonToQuark = 1', # aMC@NLO-PY8 interface settings
	    'TimeShower:limitPTmaxGlobal = on',
	    'TimeShower:nMaxGlobalRecoil = 1',
	    'TimeShower:globalRecoilMode = 2',
	    'TimeShower:nMaxGlobalBranch = 1',
            '15:offIfAny = 11 13'),
        parameterSets = cms.vstring(
            'pythia8CommonSettings', 
            'pythia8CUEP8M1Settings', 
            'processParameters')
    )
)

process.bbh_mumu_analysis = cms.EDAnalyzer("bbh_mumu",
# mH=300 GeV
# tau tau
#  	    HistOutFile = cms.untracked.string('bbh_tautau_700gev_mg5_yb2_alphadef025_pdfrwt.root'),
#	    HistOutFile = cms.untracked.string('CPtestMG5_CPmix.root'),
#	    HistOutFile = cms.untracked.string('CPtestMG5_H.root'),
#	    HistOutFile = cms.untracked.string('CPtestMG5_A.root'),
#	    HistOutFile = cms.untracked.string('bbh_tautau_700gev_mg5_yb2_alphadef025_div_sqrt2.root'),
#
#	    HistOutFile = cms.untracked.string('bbh_mumu_30gev_mg5_yb2_alpha05xdef025.root'),
#	    HistOutFile = cms.untracked.string('bbh_mumu_30gev_mg5_ybyt_alpha05xdef025.root'),
#
#	    HistOutFile = cms.untracked.string('bbh_mumu_30gev_mg5_yb2_alpha2xdef025.root'),
#	    HistOutFile = cms.untracked.string('bbh_mumu_30gev_mg5_ybyt_alpha2xdef025.root'),
#
#	    HistOutFile = cms.untracked.string('bbh_mumu_30gev_mg5_ybyt_alphadef025.root'),
#	    HistOutFile = cms.untracked.string('bbh_mumu_30gev_mg5_yb2_alphadef025_pdf.root'),
#	    HistOutFile = cms.untracked.string('bbh_mumu_30gev_mg5_ybyt_alphadef025_pdf.root'),
# Marius's events
	    HistOutFile = cms.untracked.string('bbh_mumu_30gev_mg5_yb2_alphadef025_marius.root'),
#	    HistOutFile = cms.untracked.string('bbh_mumu_30gev_mg5_ybyt_alphadef025_marius.root'),
# POWHEG SM gg->h with t+b
#	    HistOutFile = cms.untracked.string('test.root'),
#	    HistOutFile = cms.untracked.string('sm_gg_h_mumu_30gev_powheg.root'),
#           HistOutFile = cms.untracked.string('mssm_gg_a_mumu_ma30tanb40_powheg_hfact7gev.root'),
#           HistOutFile = cms.untracked.string('mssm_gg_a_mumu_ma30tanb40_powheg_hfact15gev.root'),
#	    HistOutFile = cms.untracked.string('mssm_gg_a_mumu_ma30tanb40_powheg_hfact30gev.root'),
#           HistOutFile = cms.untracked.string('mg5-amcnlo-SM_ggA_with_b_loop_only_marius_hfact15gev.root'),
#           HistOutFile = cms.untracked.string('mg5-amcnlo-SM_ggA_b_loop_g_bb_oliver-mattelaer.root'),
	    gen_jets_ak4 = cms.InputTag("ak4Jets")
)

process.load("Configuration.StandardSequences.Generator_cff")

process.load("Configuration.StandardSequences.VtxSmearedNoSmear_cff")

process.genParticles.abortOnUnknownPDGCode = False

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printList = cms.EDAnalyzer("ParticleListDrawer",
	src = cms.InputTag("genParticles"),
	maxEventsToPrint = cms.untracked.int32(-1)
)

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
	src = cms.InputTag("genParticles"),
	printP4 = cms.untracked.bool(False),
	printPtEtaPhi = cms.untracked.bool(False),
	printVertex = cms.untracked.bool(True),
	printStatus = cms.untracked.bool(False),
	printIndex = cms.untracked.bool(False),
	status = cms.untracked.vint32(1, 2, 3)
)

process.p = cms.Path(
	process.printList *
	process.printTree
)

process.load("Configuration.EventContent.EventContent_cff")

process.GEN = cms.OutputModule("PoolOutputModule",
	process.FEVTSIMEventContent,
	dataset = cms.untracked.PSet(dataTier = cms.untracked.string('GEN')),
	SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p0')),
	fileName = cms.untracked.string('test.root')
)

process.p0 = cms.Path(
	process.generator *
	process.pgen
)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.dumpevent = cms.Path(process.dump)

process.genjets04 = cms.Path(process.genParticlesForAK4Jets*process.ak4Jets)

process.an = cms.Path(process.bbh_mumu_analysis)

process.schedule = cms.Schedule(process.p0, process.genjets04, process.an)
