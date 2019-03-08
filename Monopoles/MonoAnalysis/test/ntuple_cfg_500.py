import FWCore.ParameterSet.Config as cms

process = cms.Process("Mpl")

### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")

### Standard Configurations
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
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## Fitter-smoother: loosen outlier rejection as for first data-taking with LHC "collisions"
process.KFFittingSmootherWithOutliersRejectionAndRK.BreakTrajWith2ConsecutiveMissing = False
process.KFFittingSmootherWithOutliersRejectionAndRK.EstimateCut = 1000

### Conditions
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')

### Track refitter specific stuff
#from RecoTracker.TrackProducer.TrackRefitters_cff import *
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")

### unclean EE
process.uncleanEERecovered = cms.EDProducer('UncleanSCRecoveryProducer',

            # input collections:
            cleanBcCollection = cms.InputTag('multi5x5SuperClusters','multi5x5EndcapBasicClusters'),
            cleanScCollection = cms.InputTag('multi5x5SuperClusters','multi5x5EndcapSuperClusters'),
                                    
            uncleanBcCollection = cms.InputTag('multi5x5SuperClusters','uncleanOnlyMulti5x5EndcapBasicClusters'),
            uncleanScCollection = cms.InputTag('multi5x5SuperClusters','uncleanOnlyMulti5x5EndcapSuperClusters'),
            # names of collections to be produced:
            bcCollection = cms.string('uncleanEndcapBasicClusters'),
            scCollection = cms.string('uncleanEndcapSuperClusters'),

            )

process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         #'file:/user/melsawy/Monopole_2018/CMSSW_8_0_28_patch1/src/Monopoles/MonoAnalysis/test/Work/Reco_MC_13TeV/step2_PU2016_FEVTSIM_500_100.root'


'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_100.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_101.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_102.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_103.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_104.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_105.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_106.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_107.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_108.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_109.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_110.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_111.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_112.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_113.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_114.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_115.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_116.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_117.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_118.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_119.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_120.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_121.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_122.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_123.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_124.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_125.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_126.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_127.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_128.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_129.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_130.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_131.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_132.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_133.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_134.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_135.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_136.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_137.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_138.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_139.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_140.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_141.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_142.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_143.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_144.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_145.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_146.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_147.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_148.root'
,'root://eoscms//eos/cms/store/user/srimanob/monopole/13TeV/v4/monopole_spinhalf_500GeV/step2_PU2016_FEVTSIM_500_149.root'

)

,
    duplicateCheckMode = cms.untracked.string('checkEachRealDataFile') 
)

##-------- Electron events of interest --------
process.HLTEle =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring("HLT_Photon175_v*"),     
#HLTPaths = cms.vstring("HLT_Photon175_v*", "HLT_MET200_v*","HLT_Photon165_HE10_v*", "HLT_PFMET170_HBHE_BeamHaloCleaned_v*"),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
 )


### Construct combined (clean and uncleanOnly Ecal clusters)
process.load("RecoEcal.EgammaClusterProducers.uncleanSCRecovery_cfi")

import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process.Monopoler = cms.EDAnalyzer('MonoNtupleDumper'
  ,isData = cms.bool(False)
  ,Output = cms.string("500GeV_Ntuple-ALL_TEST.root")
#  ,bits = cms.InputTag("TriggerResults","","HLT")
  ,TriggerResults = cms.InputTag("TriggerResults","","HLT")
  ,GeneratorTag = cms.InputTag("generatorSmeared","")
  ,PrimaryVertices = cms.InputTag("offlinePrimaryVertices","")                                 
  ,EcalEBRecHits = cms.InputTag("ecalRecHit","EcalRecHitsEB") 
  ,EcalEERecHits = cms.InputTag("ecalRecHit","EcalRecHitsEE") 
  ,HBHERecHits = cms.InputTag("hbhereco","")
  ,JetTag = cms.InputTag("ak4PFJets","")
  ,ElectronTag = cms.InputTag("gedGsfElectrons","")
  ,PhotonTag = cms.InputTag("photons","")
  ,METTag = cms.InputTag("pfMet","")
  ,bcClusterTag = cms.InputTag("hybridSuperClusters","uncleanOnlyHybridBarrelBasicClusters") 
  ,ccClusterTag = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters")
  ,combClusterTag = cms.InputTag("uncleanSCRecovered","uncleanHybridBarrelBasicClusters") 
  ,eeCleanTag = cms.InputTag("multi5x5SuperClusters","multi5x5EndcapBasicClusters")
  ,eeUncleanTag = cms.InputTag("multi5x5SuperClusters","uncleanOnlyMulti5x5EndcapBasicClusters") 
  ,eeCombTag = cms.InputTag("uncleanEERecovered","uncleanEndcapBasicClusters")
  ,StripSeedLength = cms.uint32(3)
  ,ClusterLength = cms.uint32(5)
  ,SeedThreshold = cms.double(50.)
  ,TrackTag=cms.InputTag("TrackRefitter")                                 
  ,TrackSource=cms.string("TrackRefitter")
  ,TrackChi2Cut=cms.untracked.double(7.5)
  ,TrackPtCut=cms.untracked.double(3.0)
  ,TrackDeDxCut=cms.untracked.double(0)
  ,TrackDefaultError=cms.untracked.double(0.05)
  ,TrackErrorFudge=cms.untracked.double(0.02)
  ,TrackHitOutput=cms.untracked.bool(True)
)

process.ecalCombine_step = cms.Path(process.uncleanSCRecovered)
process.ecalCombineEE_step = cms.Path(process.uncleanEERecovered)
#process.refit_step = cms.Path(process.TrackRefitter)
process.refit_step = cms.Path(process.MeasurementTrackerEvent * process.TrackRefitter)
process.mpl_step = cms.Path(process.Monopoler)
process.HLT_step = cms.Path(process.HLTEle)

process.options = cms.untracked.PSet(     wantSummary = cms.untracked.bool(True) )


process.p1 = cms.Schedule(process.ecalCombine_step
                          ,process.ecalCombineEE_step
                          ,process.refit_step
                          ,process.mpl_step
                          ,process.HLT_step
)
#process.outpath = cms.EndPath(process.TRACKS)
