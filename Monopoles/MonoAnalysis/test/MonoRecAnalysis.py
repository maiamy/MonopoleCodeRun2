import FWCore.ParameterSet.Config as cms

process = cms.Process("L1Example")


### Standard Configurations                                                                                                                   


#process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')




import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis') 
options.register('isMC',True,options.multiplicity.singleton,options.varType.bool," whether we are running on MC or not")
options.parseArguments()


print options.inputFiles

process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring(


'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_1000GeV/step2_PU2016_FEVTSIM_1000_100.root'

),

duplicateCheckMode = cms.untracked.string('checkEachRealDataFile') 
)



process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
reportEvery = cms.untracked.int32(5000),
limit = cms.untracked.int32(10000000)
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')



#from Configuration.AlCa.autoCond import autoCond
#from Configuration.AlCa.GlobalTag import GlobalTag
#if options.isMC:
 #   process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')
#else:
 #   process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v9 ','')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')


process.TFileService = cms.Service("TFileService",
      fileName = cms.string("Rec_L1.root"),
      closeFileFast = cms.untracked.bool(True)
  )



# set the number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.l1MenuExample = cms.EDAnalyzer("MonoRecAnalysis",
 
hltProcess=cms.string("HLT")
,Output = cms.string("MonoRec_L1_2.root")
,EcalEBRecHits = cms.InputTag("ecalRecHit","EcalRecHitsEB")
 ,GeneratorTag = cms.InputTag("generatorSmeared","")

)

process.p = cms.Path(process.l1MenuExample)
