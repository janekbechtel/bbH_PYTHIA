import FWCore.ParameterSet.Config as cms

source = cms.Source("LHESource",
	fileNames = cms.untracked.vstring('file:/localscratch/MG5_aMC_v2.2.1/MG5_aMC_v2_2_1/mumuj/Events/run_04/events.lhe'),
)

generator = cms.EDFilter("Pythia8HadronizerFilter",
	eventsToPrint = cms.untracked.uint32(2),
        UseExternalGenerators = cms.untracked.bool(True),
	maxEventsToPrint = cms.untracked.int32(2),
    	pythiaPylistVerbosity = cms.untracked.int32(2),
    	filterEfficiency = cms.untracked.double(1.0),
    	pythiaHepMCVerbosity = cms.untracked.bool(False), 
        comEnergy = cms.double(13000.),
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
                        'TimeShower:nPartonsInBorn = 2',
			'Check:epTolErr = 0.01',
        		'SLHA:keepSM = on',
        		'SLHA:minMassSM = 10.',
        		'25:onMode = off', # turn OFF all tau decays
			'25:onIfMatch = 13 -13'
			),  
	       		parameterSets = cms.vstring('processParameters'))
)

