import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = "106X_mcRun3_2024_realistic_v4"


#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         # gluino 2024
         '/store/mc/Run3Summer19DR/HSCPgluino_M_1800_TuneCP5_14TeV_pythia8/AODSIM/HSCP_SIM_106X_mcRun3_2024_realistic_v4-v2/260000/E3F00FAB-DE78-7B43-AECB-A31B6BEE8BBB.root',
         '/store/mc/Run3Summer19DR/HSCPgluino_M_1800_TuneCP5_14TeV_pythia8/AODSIM/HSCP_SIM_106X_mcRun3_2024_realistic_v4-v2/260000/7C749E19-F085-E74C-8A75-21925E5FF646.root',
         '/store/mc/Run3Summer19DR/HSCPgluino_M_1800_TuneCP5_14TeV_pythia8/AODSIM/HSCP_SIM_106X_mcRun3_2024_realistic_v4-v2/260000/74997560-6AA0-E843-AB54-9F220CE24BC1.root',
         '/store/mc/Run3Summer19DR/HSCPgluino_M_1800_TuneCP5_14TeV_pythia8/AODSIM/HSCP_SIM_106X_mcRun3_2024_realistic_v4-v2/260000/21912A8F-B46A-DA43-AD7C-756C606A6565.root'
    )
)

process.stage = cms.EDAnalyzer('ntuple'
     , format_file       = cms.string('AOD')
     , primaryVertexColl  = cms.InputTag('offlinePrimaryVertices')  #->AOD
     , isotracks             = cms.InputTag("isolatedTracks")
     , tracks             = cms.InputTag("generalTracks")
     , collectionHSCP     = cms.InputTag("HSCParticleProducer")
     , isoHSCP            = cms.InputTag("HSCPIsolation","R03")
     , dedx               = cms.InputTag("dedxHitInfo")
     , MiniDedx           = cms.InputTag("isolatedTracks")
     , dEdxHitInfoPrescale = cms.InputTag("dedxHitInfo","prescale")
     , muons             = cms.InputTag("muons")
     , muonTOF            = cms.InputTag("muons","combined")
     , muonTDT            = cms.InputTag("muons","dt")
     , muonTCSC            = cms.InputTag("muons","csc")
     , printOut           = cms.untracked.int32(-1)
     , GenPart            = cms.InputTag("genParticles")
     , runOnGS            = cms.bool(False)
     , stripSimLinks      = cms.InputTag("simSiStripDigis")
     , ROUList            = cms.vstring(
                                      'TrackerHitsTIBLowTof',  # hard-scatter (prompt) collection, module label g4simHits
                                      'TrackerHitsTIBHighTof',
                                      'TrackerHitsTIDLowTof',
                                      'TrackerHitsTIDHighTof',
                                      'TrackerHitsTOBLowTof',
                                      'TrackerHitsTOBHighTof',
                                      'TrackerHitsTECLowTof',
                                      'TrackerHitsTECHighTof'
                                      )
)



process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('/opt/sbg/cms/ui3_data1/ccollard/HSCP/gluino1800_run3/test.root')
 )

process.load("SUSYBSMAnalysis.HSCP.HSCParticleProducer_cff")
#process.p = cms.Path(process.stage)

process.p = cms.Path(process.HSCParticleProducerSeq + process.stage)
#process.p = cms.Path(process.dump)


## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.default.limit = 10
