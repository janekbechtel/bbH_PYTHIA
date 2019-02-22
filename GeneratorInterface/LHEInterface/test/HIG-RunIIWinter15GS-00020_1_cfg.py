# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/HIG-RunIIWinter15GS-00020-fragment.py --fileout file:HIG-RunIIWinter15GS-00020.root --mc --eventcontent RAWSIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --conditions MCRUN2_71_V1::All --beamspot NominalCollision2015 --step GEN,SIM --magField 38T_PostLS1 --python_filename HIG-RunIIWinter15GS-00020_1_cfg.py --no_exec -n 53
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(53)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('Configuration/GenProduction/python/HIG-RunIIWinter15GS-00020-fragment.py nevts:53'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:HIG-RunIIWinter15GS-00020.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_71_V1::All', '')

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    SLHATableForPythia8 = cms.string('## Important note!\nBLOCK SPINFO\n     1   FeynHiggs\n     2   2.10.3\n     2   built on Mar 12, 2015\nBLOCK MODSEL\n         1                  0   # Model\n         2                  1   # GridPts\n         3                  0   # Content\n         4                  0   # RPV\n         5                  0   # CPV\n         6                  0   # FV\nBLOCK SMINPUTS\n         1     1.28952896E+02   # invAlfaMZ\n         2     1.16637000E-05   # GF\n         3     1.19000000E-01   # AlfasMZ\n         4     9.11876000E+01   # MZ\n         5     4.16000000E+00   # Mb\n         6     1.72500000E+02   # Mt\n         7     1.77703000E+00   # Mtau\n        11     5.10998902E-04   # Me\n        13     1.05658357E-01   # Mmu\n        21     6.00000000E-03   # Md\n        22     3.00000000E-03   # Mu\n        23     9.50000000E-02   # Ms\n        24     1.28600000E+00   # Mc\nBLOCK MINPAR\n         3     3.00000000E+01   # TB\nBLOCK EXTPAR\n         0     0.00000000E+00   # Q\n         1     9.54716519E+01   # M1\n         2     2.00000000E+02   # M2\n         3     1.50000000E+03   # M3\n        11     1.50666667E+03   # At\n        12     1.50666667E+03   # Ab\n        13     1.50666667E+03   # Atau\n        23     2.00000000E+02   # MUE\n        25     3.00000000E+01   # TB\n        26     3.00000000E+02   # MA0\n        27     3.10586282E+02   # MHp\n        31     5.00000000E+02   # MSL(1)\n        32     5.00000000E+02   # MSL(2)\n        33     1.00000000E+03   # MSL(3)\n        34     5.00000000E+02   # MSE(1)\n        35     5.00000000E+02   # MSE(2)\n        36     1.00000000E+03   # MSE(3)\n        41     1.50000000E+03   # MSQ(1)\n        42     1.50000000E+03   # MSQ(2)\n        43     1.00000000E+03   # MSQ(3)\n        44     1.50000000E+03   # MSU(1)\n        45     1.50000000E+03   # MSU(2)\n        46     1.00000000E+03   # MSU(3)\n        47     1.50000000E+03   # MSD(1)\n        48     1.50000000E+03   # MSD(2)\n        49     1.00000000E+03   # MSD(3)\nBLOCK MASS\n   1000012     4.95834286E+02   # MSf(1,1,1)\n   1000011     5.01843810E+02   # MSf(1,2,1)\n   2000011     5.02295880E+02   # MSf(2,2,1)\n   1000002     1.49902739E+03   # MSf(1,3,1)\n   2000002     1.49958945E+03   # MSf(2,3,1)\n   1000001     1.50020509E+03   # MSf(1,4,1)\n   2000001     1.50117722E+03   # MSf(2,4,1)\n   1000014     4.95834286E+02   # MSf(1,1,2)\n   1000013     5.01398885E+02   # MSf(1,2,2)\n   2000013     5.02740033E+02   # MSf(2,2,2)\n   1000004     1.49902793E+03   # MSf(1,3,2)\n   2000004     1.49959001E+03   # MSf(2,3,2)\n   1000003     1.50017147E+03   # MSf(1,4,2)\n   2000003     1.50121083E+03   # MSf(2,4,2)\n   1000016     9.97923664E+02   # MSf(1,1,3)\n   1000015     9.97040283E+02   # MSf(1,2,3)\n   2000015     1.00502007E+03   # MSf(2,2,3)\n   1000006     8.76886935E+02   # MSf(1,3,3)\n   2000006     1.13421049E+03   # MSf(2,3,3)\n   1000005     9.92457871E+02   # MSf(1,4,3)\n   2000005     1.00955661E+03   # MSf(2,4,3)\n        25     1.25971966E+02   # Mh0\n        35     3.00344303E+02   # MHH\n        36     3.00000000E+02   # MA0\n        37     3.12323193E+02   # MHp\n   1000022     8.83213889E+01   # MNeu(1)\n   1000023     1.52049632E+02   # MNeu(2)\n   1000025    -2.10470056E+02   # MNeu(3)\n   1000035     2.65570688E+02   # MNeu(4)\n   1000024     1.48692309E+02   # MCha(1)\n   1000037     2.66117031E+02   # MCha(2)\n   1000021     1.50000000E+03   # MGl\nBLOCK DMASS\n         0     1.72500000E+02   # Q\n        25     1.36402610E+00   # Delta Mh0\n        35     6.29082132E-03   # Delta MHH\n        36     0.00000000E+00   # Delta MA0\n        37     8.41295755E-02   # Delta MHp\nBLOCK NMIX\n     1   1     9.33161766E-01   # ZNeu(1,1)\n     1   2    -1.11637757E-01   # ZNeu(1,2)\n     1   3     3.09461555E-01   # ZNeu(1,3)\n     1   4    -1.44843625E-01   # ZNeu(1,4)\n     2   1    -3.14536695E-01   # ZNeu(2,1)\n     2   2    -6.93472054E-01   # ZNeu(2,2)\n     2   3     5.12603053E-01   # ZNeu(2,3)\n     2   4    -3.96738310E-01   # ZNeu(2,4)\n     3   1     9.74524166E-02   # ZNeu(3,1)\n     3   2    -1.35722538E-01   # ZNeu(3,2)\n     3   3    -6.77905232E-01   # ZNeu(3,3)\n     3   4    -7.15909852E-01   # ZNeu(3,4)\n     4   1    -1.44148578E-01   # ZNeu(4,1)\n     4   2     6.98722344E-01   # ZNeu(4,2)\n     4   3     4.26516298E-01   # ZNeu(4,3)\n     4   4    -5.55960539E-01   # ZNeu(4,4)\nBLOCK UMIX\n     1   1    -6.06292879E-01   # UCha(1,1)\n     1   2     7.95241438E-01   # UCha(1,2)\n     2   1     7.95241438E-01   # UCha(2,1)\n     2   2     6.06292879E-01   # UCha(2,2)\nBLOCK VMIX\n     1   1    -7.95241438E-01   # VCha(1,1)\n     1   2     6.06292879E-01   # VCha(1,2)\n     2   1     6.06292879E-01   # VCha(2,1)\n     2   2     7.95241438E-01   # VCha(2,2)\nBLOCK STAUMIX\n     1   1     6.96989496E-01   # USf(1,1)\n     1   2     7.17081336E-01   # USf(1,2)\n     2   1     7.17081336E-01   # USf(2,1)\n     2   2    -6.96989496E-01   # USf(2,2)\nBLOCK STOPMIX\n     1   1     7.08257287E-01   # USf(1,1)\n     1   2    -7.05954401E-01   # USf(1,2)\n     2   1     7.05954401E-01   # USf(2,1)\n     2   2     7.08257287E-01   # USf(2,2)\nBLOCK SBOTMIX\n     1   1     6.76310146E-01   # USf(1,1)\n     1   2     7.36616988E-01   # USf(1,2)\n     2   1     7.36616988E-01   # USf(2,1)\n     2   2    -6.76310146E-01   # USf(2,2)\nBLOCK ALPHA\n              -4.61501950E-02   # Alpha\nBLOCK DALPHA\n               5.66773518E-04   # Delta Alpha\nBLOCK HMIX Q= -0.99900000E+03\n         1     2.00000000E+02   # MUE\n         2     3.00000000E+01   # TB\nBLOCK MSOFT Q=  0.00000000E+00\n         1     9.54716519E+01   # M1\n         2     2.00000000E+02   # M2\n         3     1.50000000E+03   # M3\n        31     5.00000000E+02   # MSL(1)\n        32     5.00000000E+02   # MSL(2)\n        33     1.00000000E+03   # MSL(3)\n        34     5.00000000E+02   # MSE(1)\n        35     5.00000000E+02   # MSE(2)\n        36     1.00000000E+03   # MSE(3)\n        41     1.50000000E+03   # MSQ(1)\n        42     1.50000000E+03   # MSQ(2)\n        43     1.00000000E+03   # MSQ(3)\n        44     1.50000000E+03   # MSU(1)\n        45     1.50000000E+03   # MSU(2)\n        46     1.00000000E+03   # MSU(3)\n        47     1.50000000E+03   # MSD(1)\n        48     1.50000000E+03   # MSD(2)\n        49     1.00000000E+03   # MSD(3)\nBLOCK AE Q=  0.00000000E+00\n     1   1     0.00000000E+00   # Af(1,1)\n     2   2     0.00000000E+00   # Af(2,2)\n     3   3     1.50666667E+03   # Af(3,3)\nBLOCK AU Q=  0.00000000E+00\n     1   1     0.00000000E+00   # Af(1,1)\n     2   2     0.00000000E+00   # Af(2,2)\n     3   3     1.50666667E+03   # Af(3,3)\nBLOCK AD Q=  0.00000000E+00\n     1   1     0.00000000E+00   # Af(1,1)\n     2   2     0.00000000E+00   # Af(2,2)\n     3   3     1.50666667E+03   # Af(3,3)\nBLOCK YE Q=  0.00000000E+00\n     1   1     8.80994160E-05   # Yf(1,1)\n     2   2     1.82161635E-02   # Yf(2,2)\n     3   3     3.06371119E-01   # Yf(3,3)\nBLOCK YU Q=  0.00000000E+00\n     1   1     1.72406273E-05   # Yf(1,1)\n     2   2     7.39048222E-03   # Yf(2,2)\n     3   3     9.91336068E-01   # Yf(3,3)\nBLOCK YD Q=  0.00000000E+00\n     1   1     1.00363958E-03   # Yf(1,1)\n     2   2     1.58899417E-02   # Yf(2,2)\n     3   3     6.54339723E-01   # Yf(3,3)\nBLOCK VCKMIN\n         1     2.25300000E-01   # lambda\n         2     8.08000000E-01   # A\n         3     1.32000000E-01   # rhobar\n         4     3.41000000E-01   # etabar\nBLOCK MSL2 Q=  0.00000000E+00\n     1   1     2.50000000E+05   # MSL2(1,1)\n     2   2     2.50000000E+05   # MSL2(2,2)\n     3   3     1.00000000E+06   # MSL2(3,3)\nBLOCK MSE2 Q=  0.00000000E+00\n     1   1     2.50000000E+05   # MSE2(1,1)\n     2   2     2.50000000E+05   # MSE2(2,2)\n     3   3     1.00000000E+06   # MSE2(3,3)\nBLOCK MSQ2 Q=  0.00000000E+00\n     1   1     2.25000000E+06   # MSQ2(1,1)\n     2   2     2.25000000E+06   # MSQ2(2,2)\n     3   3     1.00000000E+06   # MSQ2(3,3)\nBLOCK MSU2 Q=  0.00000000E+00\n     1   1     2.25000000E+06   # MSU2(1,1)\n     2   2     2.25000000E+06   # MSU2(2,2)\n     3   3     1.00000000E+06   # MSU2(3,3)\nBLOCK MSD2 Q=  0.00000000E+00\n     1   1     2.25000000E+06   # MSD2(1,1)\n     2   2     2.25000000E+06   # MSD2(2,2)\n     3   3     1.00000000E+06   # MSD2(3,3)\nBLOCK TE Q=  0.00000000E+00\n     1   1     0.00000000E+00   # Tf(1,1)\n     2   2     0.00000000E+00   # Tf(2,2)\n     3   3     4.61599152E+02   # Tf(3,3)\nBLOCK TU Q=  0.00000000E+00\n     1   1     0.00000000E+00   # Tf(1,1)\n     2   2     0.00000000E+00   # Tf(2,2)\n     3   3     1.49361301E+03   # Tf(3,3)\nBLOCK TD Q=  0.00000000E+00\n     1   1     0.00000000E+00   # Tf(1,1)\n     2   2     0.00000000E+00   # Tf(2,2)\n     3   3     9.85871849E+02   # Tf(3,3)\nBLOCK CVHMIX\n     1   1     9.99981566E-01   # UH(1,1)\n     1   2     6.07186184E-03   # UH(1,2)\n     1   3     0.00000000E+00   # UH(1,3)\n     2   1    -6.07186184E-03   # UH(2,1)\n     2   2     9.99981566E-01   # UH(2,2)\n     2   3     0.00000000E+00   # UH(2,3)\n     3   1     0.00000000E+00   # UH(3,1)\n     3   2     0.00000000E+00   # UH(3,2)\n     3   3     1.00000000E+00   # UH(3,3)\nBLOCK PRECOBS\n         1     4.17666251E-05   # DeltaRho\n         2     8.03586389E+01   # MWMSSM\n         3     8.03562872E+01   # MWSM\n         4     2.31529979E-01   # SW2effMSSM\n         5     2.31543011E-01   # SW2effSM\n        11     3.41956748E-09   # gminus2mu\n        21     0.00000000E+00   # EDMeTh\n        22     0.00000000E+00   # EDMn\n        23     0.00000000E+00   # EDMHg\n        31     6.02640601E-04   # bsgammaMSSM\n        32     3.98954070E-04   # bsgammaSM\n        33     2.15605651E+01   # DeltaMsMSSM\n        34     2.12958950E+01   # DeltaMsSM\n        35     4.00884141E-08   # BsmumuMSSM\n        36     3.50847021E-09   # BsmumuSM\nDECAY        25     5.94351547E-03   # Gamma(h0)\n     1.60711170E-03   2        22        22   # BR(h0 -> photon photon)\n     1.08116991E-03   2        22        23   # BR(h0 -> photon Z)\n     2.03603611E-02   2        23        23   # BR(h0 -> Z Z)\n     1.65019806E-01   2       -24        24   # BR(h0 -> W W)\n     4.71041715E-02   2        21        21   # BR(h0 -> gluon gluon)\n     6.46086023E-09   2       -11        11   # BR(h0 -> Electron electron)\n     2.87391170E-04   2       -13        13   # BR(h0 -> Muon muon)\n     8.20132866E-02   2       -15        15   # BR(h0 -> Tau tau)\n     1.36340694E-07   2        -2         2   # BR(h0 -> Up up)\n     1.88850009E-02   2        -4         4   # BR(h0 -> Charm charm)\n     1.02802495E-06   2        -1         1   # BR(h0 -> Down down)\n     2.58165131E-04   2        -3         3   # BR(h0 -> Strange strange)\n     6.63382366E-01   2        -5         5   # BR(h0 -> Bottom bottom)\nDECAY        35     3.99191224E+00   # Gamma(HH)\n     1.17290990E-06   2        22        22   # BR(HH -> photon photon)\n     1.53828836E-06   2        22        23   # BR(HH -> photon Z)\n     1.07561986E-04   2        23        23   # BR(HH -> Z Z)\n     2.40907397E-04   2       -24        24   # BR(HH -> W W)\n     3.49704730E-04   2        21        21   # BR(HH -> gluon gluon)\n     1.10043254E-08   2       -11        11   # BR(HH -> Electron electron)\n     4.89621143E-04   2       -13        13   # BR(HH -> Muon muon)\n     1.37060293E-01   2       -15        15   # BR(HH -> Tau tau)\n     7.31008328E-13   2        -2         2   # BR(HH -> Up up)\n     1.01137631E-07   2        -4         4   # BR(HH -> Charm charm)\n     1.39745180E-06   2        -1         1   # BR(HH -> Down down)\n     3.50922537E-04   2        -3         3   # BR(HH -> Strange strange)\n     8.31196607E-01   2        -5         5   # BR(HH -> Bottom bottom)\n     6.39928211E-04   2  -1000024   1000024   # BR(HH -> Chargino1 chargino1)\n     1.11390100E-02   2   1000022   1000022   # BR(HH -> neutralino1 neutralino1)\n     1.37705504E-02   2   1000022   1000023   # BR(HH -> neutralino1 neutralino2)\n     3.39505975E-03   2   1000022   1000025   # BR(HH -> neutralino1 neutralino3)\n     1.25561306E-03   2        25        25   # BR(HH -> h0 h0)\nDECAY        36     4.23356126E+00   # Gamma(A0)\n    -2.39067163E-06   2        22        22   # BR(A0 -> photon photon)\n    -7.20611526E-06   2        22        23   # BR(A0 -> photon Z)\n    -2.55169610E-04   2        21        21   # BR(A0 -> gluon gluon)\n    -1.03584641E-08   2       -11        11   # BR(A0 -> Electron electron)\n    -4.60884558E-04   2       -13        13   # BR(A0 -> Muon muon)\n     1.28979677E-01   2       -15        15   # BR(A0 -> Tau tau)\n    -2.94111404E-13   2        -2         2   # BR(A0 -> Up up)\n    -4.11620798E-08   2        -4         4   # BR(A0 -> Charm charm)\n    -1.31505131E-06   2        -1         1   # BR(A0 -> Down down)\n    -3.30230496E-04   2        -3         3   # BR(A0 -> Strange strange)\n    -7.82484735E-01   2        -5         5   # BR(A0 -> Bottom bottom)\n    -3.17586019E-02   2  -1000024   1000024   # BR(A0 -> Chargino1 chargino1)\n    -1.72897199E-02   2   1000022   1000022   # BR(A0 -> neutralino1 neutralino1)\n    -3.82970129E-02   2   1000022   1000023   # BR(A0 -> neutralino1 neutralino2)\n    -2.21551578E-05   2   1000022   1000025   # BR(A0 -> neutralino1 neutralino3)\n    -1.10849557E-04   2        23        25   # BR(A0 -> Z h0)\n    -7.24219096E-37   2        25        25   # BR(A0 -> h0 h0)\nDECAY        37     2.42100923E+00   # Gamma(Hp)\n     1.92158901E-08   2       -11        12   # BR(Hp -> Electron nu_e)\n     8.21538988E-04   2       -13        14   # BR(Hp -> Muon nu_mu)\n     2.32371024E-01   2       -15        16   # BR(Hp -> Tau nu_tau)\n     2.24060338E-06   2        -1         2   # BR(Hp -> Down up)\n     2.54164136E-05   2        -3         2   # BR(Hp -> Strange up)\n     1.43241928E-05   2        -5         2   # BR(Hp -> Bottom up)\n     1.05789837E-07   2        -1         4   # BR(Hp -> Down charm)\n     5.60764016E-04   2        -3         4   # BR(Hp -> Strange charm)\n     2.00583219E-03   2        -5         4   # BR(Hp -> Bottom charm)\n     2.43729117E-07   2        -1         6   # BR(Hp -> Down top)\n     5.70287442E-06   2        -3         6   # BR(Hp -> Strange top)\n     6.80805297E-01   2        -5         6   # BR(Hp -> Bottom top)\n     7.84098174E-02   2   1000022   1000024   # BR(Hp -> neutralino1 chargino1)\n     4.69911190E-03   2   1000023   1000024   # BR(Hp -> neutralino2 chargino1)\n     2.78276665E-04   2        24        25   # BR(Hp -> W h0)\n     1.32506646E-07   2        24        35   # BR(Hp -> W HH)\n     1.52557132E-07   2        24        36   # BR(Hp -> W A0)\nDECAY         6     1.35195755E+00   # Gamma(top)\n     1.00000000E+00   2         5        24   # BR(top -> bottom W)\n'),
    comEnergy = cms.double(13000.0),
    crossSection = cms.untracked.double(518.3),
    maxEventsToPrint = cms.untracked.int32(1),
    PythiaParameters = cms.PSet(
        pythia8CommonSettings = cms.vstring('Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = on', 
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on'),
        pythia8CUEP8M1Settings = cms.vstring('Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:pT0Ref=2.4024', 
            'MultipartonInteractions:ecmPow=0.25208', 
            'MultipartonInteractions:expPow=1.6'),
        processParameters = cms.vstring('Higgs:useBSM = on', 
            'HiggsBSM:gg2A3bbbar = on', 
            'SLHA:allowUserOverride = off', 
            'SLHA:minMassSM = 100.', 
            'PhaseSpace:mHatMin = 210.0', 
            'PhaseSpace:mHatMax = 390.0'),
        parameterSets = cms.vstring('pythia8CommonSettings', 
            'pythia8CUEP8M1Settings', 
            'processParameters')
    )
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions
