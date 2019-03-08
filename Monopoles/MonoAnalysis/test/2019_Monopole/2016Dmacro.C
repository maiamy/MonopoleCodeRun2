{
  TChain ch("monopoles");

  ch.Add("srm://maite.iihe.ac.be:8443/pnfs/iihe/cms/store/user/melsawy/SinglePhoton/CRAB3_SinglePh2016E_HLT_2019/190220_162210/0000/2016E_HLT*");


  //  ch.Add("/pnfs/iihe/cms/store/user/melsawy/SinglePhoton/CRAB3_SinglePh2016D_v1/180707_104624/0001/Data_2016D_hlt*");
  
//  ch.Add("/pnfs/iihe/cms/store/user/melsawy/SinglePhoton/CRAB3_SinglePh2016D_v1/180707_104624/0002/Data_2016D_hlt*");
  //  ch.Add("/pnfs/iihe/cms/store/user/melsawy/SinglePhoton/CRAB3_SinglePh2016D_v1/180707_104624/0003/Data_2016D_hlt*");
  //ch.Add("/pnfs/iihe/cms/store/user/melsawy/SinglePhoton/CRAB3_SinglePh2016D_v1/180707_104624/0004/Data_2016D_hlt*");
  // ch.Add("/pnfs/iihe/cms/store/user/melsawy/SinglePhoton/CRAB3_SinglePh2016D_v1/180707_104624/0005/Data_2016D_hlt*");



  gSystem->Load("Trim_C.so");

  Trim(monopoles,"Trim_2016E_00.root");

}
