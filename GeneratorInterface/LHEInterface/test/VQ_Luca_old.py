#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")

process.source = cms.Source("LHESource",
fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/VLQ_Luca_MG5/EWK_BX_8TeV__bp_MQ160GeV_WQ20.1GeV.lhe')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

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

#process.load("RecoJets.Configuration.RecoGenJets_cff")
#process.ak4PartonJets  =  process.ak4GenJets.clone()
#process.ak4PartonJets.src = cms.InputTag("genParticlesForPartonJets")

process.load("RecoJets.JetProducers.ak5GenJets_cfi")
process.ak5PartonJets  =  process.ak5GenJets.clone()
process.ak5PartonJets.src = cms.InputTag("genParticlesForPartonJets")

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.genParticlesForPartonJets = process.genParticlesForJets.clone()
#process.genParticlesForPartonJets.partonicFinalState = True
process.genParticlesForPartonJets.excludeFromResonancePids = cms.vuint32(11,12,13,14,15,16,7000001,7000021,6000007)
process.genParticlesForPartonJets.ignoreParticleIDs += cms.vuint32(11,12,13,14,15,16,7000001,7000021,6000007)

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

process.vlq_analysis = cms.EDAnalyzer("vlq",
	    HistOutFile = cms.untracked.string('EWK_BX_8TeV__bp_MQ160GeV_WQ20.1GeV.root'),
#	    parton_jets = cms.InputTag("ak4PartonJets")
	    parton_or_gen_jets = cms.InputTag("ak5PartonJets")
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

process.partonjets = cms.Path(process.genParticlesForPartonJets*process.ak5PartonJets)
#process.partonjets = cms.Path(process.genParticlesForPartonJets*process.ak4PartonJets)

process.an = cms.Path(process.vlq_analysis)

#process.schedule = cms.Schedule(process.p0, process.partonjets, process.dumpevent, process.an)

process.schedule = cms.Schedule(process.p0, process.partonjets, process.an)
