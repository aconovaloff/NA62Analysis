#include "Parameters.hh"
#include "AnalysisTools.hh"
#include "Event.hh"



Parameters::Parameters(Int_t run){
if(run == 3809){
         tCedarLow =-1;//
         tCedarHi =0.5;//
         tCedarmid =-0.13; //
//         tLAVLow = -2;
//         tLAVHi = 2;
         tLAVLow = -5; //loosened LAV cuts to see if efficiency improves
         tLAVHi = 5;
//         tMUV3Low = -2;
//         tMUV3Hi =2;
 tMUV3Low = -5;
         tMUV3Hi =5;
//         tCHANTILow = -4;
//         tCHANTIHi = 2;
         tCHANTILow = -5;
         tCHANTIHi = 5;
         tMUV1Low = -10;
         tMUV1Hi = 10;
         tMUV2Low = -5;
         tMUV2Hi = 5;
         tRICHLow = -1;
	 tRICHHi = 1;
//         tRICHHi = 0.5;
//	 tLKrLow = -4;
//	 tLKrHi = 2;
         tLKrLow = -5;
         tLKrHi = 5;
//	 tLKrLow = -15;//per G's recommendation
  //       tLKrHi = 10;//per G's recommendation
         tStrawLow = -15;
         tStrawHi = 15; 
         tSACLow = -10;
         tSACHi = 5;
         tSAVLow = -15;
         tSAVHi = 5;
         tIRCLow = -8;
         tIRCHi = 5;
	 tGTKLow = -4;
	 tGTKHi = 4;
}//////////////these use CHOD as reference time
         LKrStartPos = 240413;
         MUV3StartPos = 246850;
         RICHStartPos = 236875;
         MUV1StartPos = 244341;
	MUV2StartPos = 245290;
         CHODStartPos = 239009;
         mK = 0.493677; //in GeV
         mPi = 0.13957018; //in GeV


}

