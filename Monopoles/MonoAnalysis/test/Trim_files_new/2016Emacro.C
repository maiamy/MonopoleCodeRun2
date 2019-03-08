{
  TChain ch("monopoles");
  ch.Add("/pnfs/iihe/cms/store/user/melsawy/SinglePhoton/CRAB3_SinglePh2016E_HLT_2019/190220_162210/0000/2016E_HLT*");

  gSystem->Load("Trim_C.so");
  Trim(monopoles,"Trim_E_00.root");
}
