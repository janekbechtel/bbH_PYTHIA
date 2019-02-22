#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")

process.source = cms.Source("LHESource",
#	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_6_1/CPggH/Events /run_17/events.lhe')
#	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_6_1/CPggH/Events/run_18_LO/events.lhe')
	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_6_1/CPggHLO/Events/run_01/unweighted_events.lhe')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.options = cms.untracked.PSet(
SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.configurationMetadata = cms.untracked.PSet(
	version = cms.untracked.string('alpha'),
	name = cms.untracked.string('LHEF input'),
	annotation = cms.untracked.string('h2j'),
)

process.load("Configuration.StandardSequences.Services_cff")

process.RandomNumberGeneratorService.generator = cms.PSet(
	initialSeed = cms.untracked.uint32(123456789),
	engineName = cms.untracked.string('HepJamesRandom')
)

from Configuration.Generator.Pythia8CommonSettings_cfi import *

from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

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
              pythia8CommonSettingsBlock,
              pythia8CUEP8M1SettingsBlock,
              processParameters = cms.vstring(
                'JetMatching:setMad = off',
            	'JetMatching:scheme = 1',
            	'JetMatching:merge = on',
            	'JetMatching:jetAlgorithm = 2',
            	'JetMatching:etaJetMax = 5.',
            	'JetMatching:coneRadius = 1.',
            	'JetMatching:slowJetPower = 1',
            	'JetMatching:qCut = 30.', #this is the actual merging scale
            	'JetMatching:nQmatch = 5', #4 corresponds to 4-flavour scheme (no matching of b-quarks)
            	'JetMatching:nJetMax = 2', #number of partons in born matrix element for highest multiplicity
            	'JetMatching:doShowerKt = off', #off for MLM matching, turn on for shower-kT matching
        	),
              parameterSets = cms.vstring('pythia8CommonSettings',
					  'pythia8CUEP8M1Settings',
				    	  'processParameters',
				    )
    )
)

process.cpggh_analysis = cms.EDAnalyzer("cpggh",
#	    HistOutFile = cms.untracked.string('CPtestMG5_CPmix.root'),
#	    HistOutFile = cms.untracked.string('CPtestMG5_H.root'),
	    HistOutFile = cms.untracked.string('CPtestMG5_H_LO.root'),
#	    HistOutFile = cms.untracked.string('CPtestMG5_A.root'),
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

process.an = cms.Path(process.cpggh_analysis)

process.schedule = cms.Schedule(process.p0, process.genjets04, process.an)
