#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")

process.source = cms.Source("LHESource",
	  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/public/MG5_aMC_v2_6_1/CPMaxMixggH2j/Events/run_05/events.lhe')
)

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))

process.configurationMetadata = cms.untracked.PSet(
	version = cms.untracked.string('alpha'),
	name = cms.untracked.string('LHEF input'),
	annotation = cms.untracked.string('h2jets')
)

process.load("Configuration.StandardSequences.Services_cff")

process.RandomNumberGeneratorService.generator = cms.PSet(
	initialSeed = cms.untracked.uint32(123456789),
	engineName = cms.untracked.string('HepJamesRandom')
)


from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

#from Configuration.Generator.Pythia8CommonSettings_cfi import *
#from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

pythia8CommonSettingsBlock = cms.PSet(
    pythia8CommonSettings = cms.vstring(
      'Tune:preferLHAPDF = 2',
      'Main:timesAllowErrors = 10000',
      'Check:epTolErr = 0.01',
      'Beams:setProductionScalesFromLHEF = off',
      'SLHA:keepSM = on',
      'SLHA:minMassSM = 1000.',
      'ParticleDecays:limitTau0 = on',
      'ParticleDecays:tau0Max = 10',
      'ParticleDecays:allowPhotonRadiation = on',
    )
)

pythia8CP5SettingsBlock = cms.PSet(
    pythia8CP5Settings = cms.vstring(
        'Tune:pp 14',
        'Tune:ee 7',
        'MultipartonInteractions:ecmPow=0.03344',
        'PDF:pSet=20',
        'MultipartonInteractions:bProfile=2',
        'MultipartonInteractions:pT0Ref=1.41',
        'MultipartonInteractions:coreRadius=0.7634',
        'MultipartonInteractions:coreFraction=0.63',
        'ColourReconnection:range=5.176',
        'SigmaTotal:zeroAXB=off',
        'SpaceShower:alphaSorder=2',
        'SpaceShower:alphaSvalue=0.118',
        'SigmaProcess:alphaSvalue=0.118',
        'SigmaProcess:alphaSorder=2',
        'MultipartonInteractions:alphaSvalue=0.118',
        'MultipartonInteractions:alphaSorder=2',
        'TimeShower:alphaSorder=2',
        'TimeShower:alphaSvalue=0.118',
        )
)

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(2),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    UseExternalGenerators = cms.untracked.bool(True),
    ExternalDecays = cms.PSet(
        Tauola = cms.untracked.PSet(
            UseTauolaPolarization = cms.bool(True),
            InputCards = cms.PSet(
                mdtau = cms.int32(0),
                pjak2 = cms.int32(4),
                pjak1 = cms.int32(4)
#       call msg_message ("   JAK =  0    : All decay mode")
#       call msg_message ("   JAK =  1    : electron mode")
#       call msg_message ("   JAK =  2    : muon mode")
#       call msg_message ("   JAK =  3    : pion mode")
#       call msg_message ("   JAK =  4    : rho mode")
#       call msg_message ("   JAK =  5    : a1 mode")
#       call msg_message ("   JAK =  6    : K mode")
#       call msg_message ("   JAK =  7    : K* mode")
#       call msg_message ("   JAK =  8-13 : n pion modes")
#       call msg_message ("   JAK = 14-19 : K K pi and K pi pi modes")
#       call msg_message ("   JAK = 20-21 : eta pi pi; gamma pi pi modes")
            ),
            parameterSets = cms.vstring('setHiggsScalarPseudoscalarPDG', 
	                                 'setHiggsScalarPseudoscalarMixingAngle'),
            setHiggsScalarPseudoscalarPDG = cms.int32(25),
            setHiggsScalarPseudoscalarMixingAngle = cms.double(90), #this is for SM decay
                                                                   #change to 90 for pure pseudoscalar
        ),
        parameterSets = cms.vstring('Tauola'),
    ),
    PythiaParameters = cms.PSet(
        pythia8CommonSettings = cms.vstring(
    	  'Tune:preferLHAPDF = 2',
      	  'Main:timesAllowErrors = 10000',
      	  'Check:epTolErr = 0.01',
      	  'Beams:setProductionScalesFromLHEF = off',
      	  'SLHA:keepSM = on',
      	  'SLHA:minMassSM = 10.',
      	  'ParticleDecays:limitTau0 = on',
      	  'ParticleDecays:tau0Max = 10',
      	  'ParticleDecays:allowPhotonRadiation = on'),
     	pythia8CP5Settings = cms.vstring(
           'Tune:ee 7',
           'MultipartonInteractions:ecmPow=0.03344',
#           'PDF:pSet=20',
           'MultipartonInteractions:bProfile=2',
           'MultipartonInteractions:pT0Ref=1.41',
           'MultipartonInteractions:coreRadius=0.7634',
           'MultipartonInteractions:coreFraction=0.63',
           'ColourReconnection:range=5.176',
           'SigmaTotal:zeroAXB=off',
           'SpaceShower:alphaSorder=2',
           'SpaceShower:alphaSvalue=0.118',
           'SigmaProcess:alphaSvalue=0.118',
           'SigmaProcess:alphaSorder=2',
           'MultipartonInteractions:alphaSvalue=0.118',
           'MultipartonInteractions:alphaSorder=2',
           'TimeShower:alphaSorder=2',
           'TimeShower:alphaSvalue=0.118'),
        processParameters = cms.vstring(
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
      		'TimeShower:weightGluonToQuark = 1',
            '25:m0 = 125.0',
            '25:onMode = off', 
            '25:onIfMatch = 15 -15'),
#            '25:onIfMatch = 13 -13'),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'processParameters')
  )
)

process.cpggh_analysis = cms.EDAnalyzer("cpggh_tautau",
	    HistOutFile = cms.untracked.string('cpoddgghntpl.root')
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

process.p0 = cms.Path(process.generator*process.pgen
)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.dumpevent = cms.Path(process.dump)

process.an = cms.Path(process.cpggh_analysis)

process.schedule = cms.Schedule(process.p0, process.an)
