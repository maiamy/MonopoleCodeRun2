from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

my_cern_username = getUsernameFromSiteDB()
config.section_('General')
config.General.requestName = 'SP2016B_HLT_v1_3'
config.General.workArea = 'DataSP2016B_v1_3'
config.General.transferOutputs = True
#config.General.transferLogs = True

config.section_('JobType')
config.JobType.psetName = 'ntuple2016B_v1.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['2016B_HLT_v1.root']

config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/group/phys_exotica/Monopole2017'
config.Data.publishWithGroupName = False
config.Data.publication = False
config.Data.inputDataset ='/SinglePhoton/Run2016B-EXOMONOPOLE-07Aug17_ver1-v1/USER'
config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.splitting = 'Automatic'
config.Data.outputDatasetTag = 'CRAB3_SinglePh2016B_2019_v1'


#config.section_('User')
config.section_('Site')

config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_BE_IIHE' 
