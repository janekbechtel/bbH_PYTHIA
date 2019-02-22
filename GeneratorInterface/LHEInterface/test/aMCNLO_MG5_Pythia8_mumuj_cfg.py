# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: myfragment.py -s GEN --no_exec --conditions auto:mc --eventcontent RAWSIM
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("LHESource",
#  fileNames = cms.untracked.vstring('file:/localscratch/MG5_aMC_v2.2.1/MG5_aMC_v2_2_1/nunujj_nlo/Events/run_07/events.lhe')
#  fileNames = cms.untracked.vstring('file:/localscratch/MG5_aMC_v2.2.1/MG5_aMC_v2_2_1/mumujj_nlo/Events/run_06/events.lhe')
  fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/anikiten/localscratch/MG5_aMC_v2.2.1/MG5_aMC_v2_2_1/mumuj/Events/run_05/events.lhe')
#  fileNames = cms.untracked.vstring('file:/localscratch/MG5_aMC_v2_2_3/mumubb_nlo/Events/run_04/events.lhe')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('myfragment.py nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('myfragment_py_GEN.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    UseExternalGenerators = cms.untracked.bool(True),
    pythiaPylistVerbosity = cms.untracked.int32(2),
    filterEfficiency = cms.untracked.double(1.0),
    eventsToPrint = cms.untracked.uint32(2),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.0),
    maxEventsToPrint = cms.untracked.int32(2),
    PythiaParameters = cms.PSet(
        processParameters = cms.vstring('Main:timesAllowErrors = 10000', 
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
            '25:onMode = off', 
            '25:onIfMatch = 13 -13'),
        parameterSets = cms.vstring('processParameters')
    )
)

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.genParticlesForPartonJets = process.genParticlesForJets.clone()
process.genParticlesForPartonJets.partonicFinalState = True
process.genParticlesForPartonJets.excludeFromResonancePids = cms.vuint32(11, 12, 13, 14, 15, 16)

process.load("RecoJets.JetProducers.ak5GenJets_cfi")
process.ak5PartonJets  =  process.ak5GenJets.clone()
process.ak5PartonJets.src = cms.InputTag("genParticlesForPartonJets")

process.bbh_mumu_analysis = cms.EDAnalyzer("mumuj_mg5",
#     HistOutFile = cms.untracked.string('nunujj_nlo_mg5.root'),
#     HistOutFile = cms.untracked.string('mumujj_nlo_mg5.root'),
      HistOutFile = cms.untracked.string('test.root'),
##     HistOutFile = cms.untracked.string('mumuj_mg5.root'),
#    HistOutFile = cms.untracked.string('mumubb_mg5_5M.root'),
    parton_jets = cms.InputTag("ak5PartonJets")
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# my staff

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.mydump = cms.Path(process.dump)

process.partonjets = cms.Path(process.genParticlesForPartonJets*process.ak5PartonJets)

process.an = cms.Path(process.bbh_mumu_analysis)

process.schedule = cms.Schedule(process.generation_step, process.partonjets, process.an,process.genfiltersummary_step, process.endjob_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 
