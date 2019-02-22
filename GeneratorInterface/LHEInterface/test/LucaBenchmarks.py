#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")

process.source = cms.Source("LHESource",
#fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/VLQ_Luca_MG5/KenLane13TeV_Zj.lhe')
#fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/VLQ_Luca_MG5/EWK_BX_8TeV__bp_MQ160GeV_WQ20.1GeV.lhe')
fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/VLQ_Luca_MG5/EWK_BX_13TeV__bp_MQ160GeV_WQ20.1GeV.lhe')
)


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))

# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000000))

# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(50000))

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

# gen jets ak4 and ak5
#
process.load("RecoJets.JetProducers.ak4GenJets_cfi")
process.ak4Jets  =  process.ak4GenJets.clone()
process.ak4Jets.src = cms.InputTag("genParticlesForAK4Jets")

process.load("RecoJets.JetProducers.ak5GenJets_cfi")
process.ak5Jets  =  process.ak5GenJets.clone()
process.ak5Jets.src = cms.InputTag("genParticlesForAK5Jets")

#input particles for gen jets 
#
process.load("RecoJets.Configuration.GenJetParticles_cff")

process.genParticlesForAK4Jets = process.genParticlesForJets.clone()
process.genParticlesForAK4Jets.excludeFromResonancePids = cms.vuint32(11,12,13,14,15,16,7000021,6000007)
process.genParticlesForAK4Jets.ignoreParticleIDs +=       cms.vuint32(11,12,13,14,15,16,7000021,6000007)

process.genParticlesForAK5Jets = process.genParticlesForJets.clone()
process.genParticlesForAK5Jets.excludeFromResonancePids = cms.vuint32(11,12,13,14,15,16,7000021,6000007)
process.genParticlesForAK5Jets.ignoreParticleIDs +=       cms.vuint32(11,12,13,14,15,16,7000021,6000007)

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
	eventsToPrint = cms.untracked.uint32(2),
        UseExternalGenerators = cms.untracked.bool(True),
	maxEventsToPrint = cms.untracked.int32(2),
    	pythiaPylistVerbosity = cms.untracked.int32(1),
    	filterEfficiency = cms.untracked.double(1.0),
    	pythiaHepMCVerbosity = cms.untracked.bool(False), 
        comEnergy = cms.double(13000.),
#	SLHAFileForPythia8 = cms.string('powheg-fh-output_tautau.slha'),
       	PythiaParameters = cms.PSet(
                processParameters = cms.vstring(
	        	'Main:timesAllowErrors = 10000',
        	        'ParticleDecays:limitTau0 = on',
			'ParticleDecays:tau0Max = 10',
        		'Tune:ee 3',
        		'Tune:pp 5',
			'SpaceShower:pTmaxMatch = 1',
			'SpaceShower:pTmaxFudge = 1',
			'SpaceShower:MEcorrections = off',
			'TimeShower:pTmaxMatch = 1',
			'TimeShower:pTmaxFudge = 1',
			'TimeShower:MEcorrections = off',
			'TimeShower:globalRecoil = on',
			'TimeShower:limitPTmaxGlobal = on',
			'TimeShower:nMaxGlobalRecoil = 1',
			'TimeShower:globalRecoilMode = 2',
			'TimeShower:nMaxGlobalBranch = 1',
			'Check:epTolErr = 0.01',
        		'SLHA:keepSM = on',
        		'SLHA:minMassSM = 10.',
#        		'25:onMode = off',
#			'25:onIfMatch = 15 -15',
#        		'36:onMode = off',
#      			'36:onIfMatch = 15 -15',
                        '15:offIfAny = 11 13'
			),  
	       		parameterSets = cms.vstring('processParameters'))
)

process.vlq_analysis = cms.EDAnalyzer("cone04vs05",
#            HistOutFile = cms.untracked.string('EWK_BX_8TeV__bp_MQ160GeV_WQ20.1GeV.root'),
	    HistOutFile = cms.untracked.string('EWK_BX_13TeV__bp_MQ160GeV_WQ20.1GeV.root'),
	    gen_jets_ak4 = cms.InputTag("ak4Jets"),
	    gen_jets_ak5 = cms.InputTag("ak5Jets")
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
process.genjets05 = cms.Path(process.genParticlesForAK5Jets*process.ak5Jets)
process.an = cms.Path(process.vlq_analysis)

process.schedule = cms.Schedule(process.p0, process.genjets04, process.genjets05, process.an)
