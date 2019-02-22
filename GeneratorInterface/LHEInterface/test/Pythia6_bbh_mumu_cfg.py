# Auto generated configuration file
# using: 
# Revision: 1.303.2.7 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: QCD_Pt_80_120_7TeV.cfi -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT:GRun,RAW2DIGI,L1Reco -n 10 --geometry DB --conditions auto:mc --relval 9000,50 --datatier GEN-SIM-DIGI-RAW-HLTDEBUG --eventcontent FEVTDEBUGHLT --no_exec --python_filename RelValQCD_Pt_80_120_MC.py
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.Generator.PythiaUEZ2Settings_cfi import *
from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *



process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(10000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('QCD_Pt_80_120_7TeV.cfi nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('QCD_Pt_80_120_7TeV_cfi_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-HLTDEBUG')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START53_V7A::All'

# parton jet
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.genParticlesForPartonJets = process.genParticlesForJets.clone()
# process.genParticlesForPartonJets.ignoreParticleIDs = cms.vuint32()
process.genParticlesForPartonJets.partonicFinalState = True
process.genParticlesForPartonJets.excludeFromResonancePids = cms.vuint32(11, 12, 13, 14, 15, 16)

process.load("RecoJets.JetProducers.ak5GenJets_cfi")
process.ak5PartonJets  =  process.ak5GenJets.clone()
process.ak5PartonJets.src = cms.InputTag("genParticlesForPartonJets")

process.nmssmanalysis = cms.EDAnalyzer("nmssm_bba2mu",
    HistOutFile = cms.untracked.string('bbh_mumu_400gev_py6.root'),
    parton_jets = cms.InputTag("ak5PartonJets")
)

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    comEnergy = cms.double(8000.0),
    crossSection = cms.untracked.double(1.1),
    ExternalDecays = cms.PSet(
        Tauola = cms.untracked.PSet(
            TauolaPolar,
            TauolaDefaultInputCards
#	     InputCards = cms.PSet(
#               pjak1 = cms.int32(0),
#               pjak2 = cms.int32(0),
#               mdtau = cms.int32(116)
#             )
        ),
        parameterSets = cms.vstring('Tauola')
    ),
    UseExternalGenerators = cms.untracked.bool(True),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring('MSEL=0         ! User defined processes',
            'MSUB(186)= 1   ! gg->QQbarA (MSSM)',
            'KFPR(186,2)= 5 ! Q = b',  
            # MSSM settings
            'IMSS(4)= 2     ! masses fixed by user',
            'RMSS(5)= 400.  ! tan beta', 
            'RMSS(19)= 30. ! m_A', 
            'RMSS(1)= 100.  ! M1',   
            'RMSS(2)= 200.  ! M2',   
            'RMSS(3)= 800.  ! Mg',   
            'RMSS(4)= 200.  ! mu',   
            'RMSS(6)= 1000.  ! MS',
            'RMSS(7)= 1000.  ! MS',
            'RMSS(8)= 1000.  ! MS',
            'RMSS(9)= 1000.  ! MS',
            'RMSS(10)= 1000.  ! MS',
            'RMSS(11)= 1000.  ! MS',
            'RMSS(12)= 1000.  ! MS',
            'RMSS(13)= 1000.  ! MS',
            'RMSS(14)= 1000.  ! MS',
            'RMSS(15)= 2000.  ! Ab',
            'RMSS(16)= 2000.  ! At',
            'RMSS(17)= 2000.  ! Atau',
            # Higgs masses
            # 'PMAS(25,1)=125.   ! mh',
            # 'PMAS(35,1)=125.  ! mH',
            'PMAS(36,1)=400.  ! mA',
            # Switch off / on desirable channels for A->tautau
            'MDME(420,1)=0  ! Higgs(H) decay into d              dbar', 
            'MDME(421,1)=0 ', 
            'MDME(422,1)=0 ', 
            'MDME(423,1)=0 ', 
            'MDME(424,1)=0 ', 
            'MDME(425,1)=0 ', 
            'MDME(426,1)=0 ', 
            'MDME(427,1)=0 ', 
            'MDME(428,1)=0 ', 
            'MDME(429,1)=1  ! decay to mu+mu-', 
            'MDME(430,1)=0  ', 
            'MDME(431,1)=0 ', 
            'MDME(432,1)=0 ', 
            'MDME(433,1)=0 ', 
            'MDME(434,1)=0 ', 
            'MDME(435,1)=0 ', 
            'MDME(436,1)=0 ', 
            'MDME(437,1)=0 ', 
            'MDME(438,1)=0 ', 
            'MDME(439,1)=0 ', 
            'MDME(440,1)=0 ', 
            'MDME(441,1)=0 ', 
            'MDME(442,1)=0 ', 
            'MDME(443,1)=0 ', 
            'MDME(444,1)=0 ', 
            'MDME(445,1)=0 ', 
            'MDME(446,1)=0 ', 
            'MDME(447,1)=0 ', 
            'MDME(448,1)=0 ', 
            'MDME(449,1)=0 ', 
            'MDME(450,1)=0 ',
            'MDME(451,1)=0 ', 
            'MDME(452,1)=0 ', 
            'MDME(453,1)=0 ', 
            'MDME(454,1)=0 ', 
            'MDME(455,1)=0 ', 
            'MDME(456,1)=0 ', 
            'MDME(457,1)=0 ', 
            'MDME(458,1)=0 ', 
            'MDME(459,1)=0 ', 
            'MDME(460,1)=0 ', 
            'MDME(461,1)=0 ', 
            'MDME(462,1)=0 ', 
            'MDME(463,1)=0 ', 
            'MDME(464,1)=0 ', 
            'MDME(465,1)=0 ', 
            'MDME(466,1)=0 ', 
            'MDME(467,1)=0 ', 
            'MDME(468,1)=0 ', 
            'MDME(469,1)=0 ', 
            'MDME(470,1)=0 ', 
            'MDME(471,1)=0 ', 
            'MDME(472,1)=0 ', 
            'MDME(473,1)=0 ', 
            'MDME(474,1)=0 ', 
            'MDME(475,1)=0 ', 
            'MDME(476,1)=0 ', 
            'MDME(477,1)=0 ', 
            'MDME(478,1)=0 ', 
            'MDME(479,1)=0 ', 
            'MDME(480,1)=0 ', 
            'MDME(481,1)=0 ', 
            'MDME(482,1)=0 ', 
            'MDME(483,1)=0 ', 
            'MDME(484,1)=0 ', 
            'MDME(485,1)=0 ', 
            'MDME(486,1)=0 ', 
            'MDME(487,1)=0 ', 
            'MDME(488,1)=0 ', 
            'MDME(489,1)=0 ', 
            'MDME(490,1)=0 ', 
            'MDME(491,1)=0 ', 
            'MDME(492,1)=0 ', 
            'MDME(493,1)=0 ', 
            'MDME(494,1)=0 ', 
            'MDME(495,1)=0 ', 
            'MDME(496,1)=0 ', 
            'MDME(497,1)=0 ', 
            'MDME(498,1)=0 ', 
            'MDME(499,1)=0 ', 
            'MDME(500,1)=0 ', 
            'MDME(501,1)=0 ', 
            'MDME(502,1)=0 '),
        # This is a vector of ParameterSet names to be read, in this order
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
# process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)

process.an = cms.Path(process.genParticlesForPartonJets*process.ak5PartonJets*process.nmssmanalysis)

#process.an = cms.Path(process.nmssmanalysis)

process.schedule = cms.Schedule(process.generation_step,process.an)

#process.p1 = cms.Path(process.genParticlesForPartonJets*process.ak5PartonJets*process.bbhanalysis)

#process.schedule.extend(process.HLTSchedule)
#process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.endjob_step,process.FEVTDEBUGHLToutput_step])
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 
