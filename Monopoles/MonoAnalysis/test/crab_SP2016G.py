from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

my_cern_username = getUsernameFromSiteDB()
config.section_('General')
config.General.requestName = 'SP201G_NEW_3'
config.General.workArea = 'DataSP2016G_HLT_new_3'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.psetName = 'ntuple2016G.py'
config.JobType.pluginName = 'Analysis'
#config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.JobType.outputFiles = ['2016G_HLT.root']

config.section_('Data')
config.Data.outLFNDirBase = '/store/group/phys_exotica/Monopole2017'
#config.Data.outLFNDirBase = '/pnfs/iihe/cms/store/user/melsawy/New_2016D'
config.Data.publishWithGroupName = False
config.Data.publication = False
#config.Data.ignoreLocality = True 
config.Data.inputDataset ='/SinglePhoton/Run2016G-EXOMONOPOLE-07Aug17-v1/USER'

config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

#config.Data.splitting = 'FileBased'
#config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3_SinglePh2016G_new'


config.section_('User')
config.section_('Site')
#config.Site.whitelist = 'T1_UK_*'
#config.Site.whitelist = ["T1_UK*"]
#config.Site.blacklist = 'T2_IT_Bari'
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_BE_IIHE' 
