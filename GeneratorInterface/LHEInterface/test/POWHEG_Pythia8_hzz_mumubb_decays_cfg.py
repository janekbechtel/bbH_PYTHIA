#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")

process.source = cms.Source("LHESource",
   	fileNames = cms.untracked.vstring('file:pwgevents_125_renf_05mh.lhe')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

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

#process.generator = cms.EDFilter("LHEProducer",
process.generator = cms.EDFilter("Pythia8HadronizerFilter",
	eventsToPrint = cms.untracked.uint32(5),

#	hadronisation = cms.PSet(
#		generator = cms.string('Pythia8'),
#		maxEventsToPrint = cms.untracked.int32(2),
#		pythiaPylistVerbosity = cms.untracked.int32(2),
#		parameterSets = cms.vstring('pythiaCMSDefaults'),
#		pythiaCMSDefaults = cms.vstring('Check:event = off')
#	),

#        ExternalDecays = cms.PSet(
#                Tauola = cms.untracked.PSet(TauolaPolar,TauolaDefaultInputCards),
#                parameterSets = cms.vstring('Tauola')),
	
        UseExternalGenerators = cms.untracked.bool(True),
	maxEventsToPrint = cms.untracked.int32(0),
    	pythiaPylistVerbosity = cms.untracked.int32(0),
    	filterEfficiency = cms.untracked.double(1.0),
    	pythiaHepMCVerbosity = cms.untracked.bool(False), 
        comEnergy = cms.double(8000.),
#        SLHAFileForPythia8 = cms.string('GeneratorInterface/LHEInterface/test/powheg-fh-output.slha'),
       	PythiaParameters = cms.PSet(
                processParameters = cms.vstring(
	        	'Main:timesAllowErrors = 10000',
        	        'ParticleDecays:limitTau0 = on',
			'ParticleDecays:tau0Max = 10',
        		'Tune:ee 3',
        		'Tune:pp 5',
        		'SLHA:keepSM = on',
        		'SLHA:minMassSM = 10.',
        		'25:onMode = off', # turn OFF all Higgs decays
			'25:onIfMatch = 23 23',
                        '23:onMode = off',
                        '23:onIfAny = 5 11 13 15',
			),  
	       		parameterSets = cms.vstring('processParameters'))
)

process.nmssmanalysis = cms.EDAnalyzer("h_zz_llbb",
    HistOutFile = cms.untracked.string('pwgevents_125_renf_05mh.root'),
)

process.load("Configuration.StandardSequences.Generator_cff")

process.p0 = cms.Path(
	process.generator *
	process.pgen
)

process.an = cms.Path(process.nmssmanalysis)

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

# process.outpath = cms.EndPath(process.GEN)

# process.schedule = cms.Schedule(process.p0, process.outpath)

process.schedule = cms.Schedule(process.p0, process.an)
