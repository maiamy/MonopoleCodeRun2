from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

my_cern_username = getUsernameFromSiteDB()
config.section_('General')
config.General.requestName = 'SP2016B_HLT'
config.General.workArea = 'DataSP2016B_v2'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.psetName = 'ntuple2016B_v2.py'
config.JobType.pluginName = 'Analysis'
#config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.JobType.outputFiles = ['2016B_HLT_v2.root']

config.section_('Data')
config.Data.outLFNDirBase = '/store/group/phys_exotica/Monopole2017'
#config.Data.outLFNDirBase = '/eos/cms/store/group/phys_exotica/Monopole_2017'
config.Data.publishWithGroupName = False
config.Data.publication = True
#config.Data.ignoreLocality = True 
config.Data.inputDataset ='/SinglePhoton/Run2016B-EXOMONOPOLE-07Aug17_ver2-v1/USER'

config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1

config.Data.splitting = 'Automatic'
#config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3_SinglePh2016B_2019_v2'


config.section_('User')
config.section_('Site')
#config.Site.whitelist = 'T1_UK_*'
#config.Site.whitelist = ["T1_UK*"]
#config.Site.blacklist = ['T1_UK_RAL_Disk']
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_BE_IIHE' 
