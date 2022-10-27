/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2017, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     TEncTop.h
    \brief    encoder class (header)
*/

#ifndef __TENCTOP__
#define __TENCTOP__

// Include files
#include "TLibCommon/TComList.h"
#include "TLibCommon/TComPrediction.h"
#include "TLibCommon/TComTrQuant.h"
#include "TLibCommon/TComLoopFilter.h"
#include "TLibCommon/AccessUnit.h"

#include "TLibVideoIO/TVideoIOYuv.h"

#include "TEncCfg.h"
#include "TEncGOP.h"
#include "TEncSlice.h"
#include "TEncEntropy.h"
#include "TEncCavlc.h"
#include "TEncSbac.h"
#include "TEncSearch.h"
#include "TEncSampleAdaptiveOffset.h"
#include "TEncPreanalyzer.h"
#include "TEncRateCtrl.h"
//! \ingroup TLibEncoder
//! \{

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// encoder class
class TEncTop : public TEncCfg
{
private:
  // picture
  Int                     m_iPOCLast;                     ///< time index (POC)
  Int                     m_iNumPicRcvd;                  ///< number of received pictures
  UInt                    m_uiNumAllPicCoded;             ///< number of coded pictures
  TComList<TComPic*>      m_cListPic;                     ///< dynamic list of pictures

  // encoder search
  TEncSearch              m_cSearch;                      ///< encoder search class
  //TEncEntropy*            m_pcEntropyCoder;                     ///< entropy encoder
  TEncCavlc*              m_pcCavlcCoder;                       ///< CAVLC encoder
  // coding tool
  TComTrQuant             m_cTrQuant;                     ///< transform & quantization class
  TComLoopFilter          m_cLoopFilter;                  ///< deblocking filter class
  TEncSampleAdaptiveOffset m_cEncSAO;                     ///< sample adaptive offset class
  TEncEntropy             m_cEntropyCoder;                ///< entropy encoder
  TEncCavlc               m_cCavlcCoder;                  ///< CAVLC encoder
  TEncSbac                m_cSbacCoder;                   ///< SBAC encoder
  TEncBinCABAC            m_cBinCoderCABAC;               ///< bin coder CABAC

  // processing unit
  TEncGOP                 m_cGOPEncoder;                  ///< GOP encoder
  TEncSlice               m_cSliceEncoder;                ///< slice encoder
  TEncCu                  m_cCuEncoder;                   ///< CU encoder
  // SPS
  ParameterSetMap<TComSPS> m_spsMap;                      ///< SPS. This is the base value. This is copied to TComPicSym
  ParameterSetMap<TComPPS> m_ppsMap;                      ///< PPS. This is the base value. This is copied to TComPicSym
  // RD cost computation
  TComRdCost              m_cRdCost;                      ///< RD cost computation class
  TEncSbac***             m_pppcRDSbacCoder;              ///< temporal storage for RD computation
  TEncSbac                m_cRDGoOnSbacCoder;             ///< going on SBAC model for RD stage
#if FAST_BIT_EST
  TEncBinCABACCounter***  m_pppcBinCoderCABAC;            ///< temporal CABAC state storage for RD computation
  TEncBinCABACCounter     m_cRDGoOnBinCoderCABAC;         ///< going on bin coder CABAC for RD stage
#else
  TEncBinCABAC***         m_pppcBinCoderCABAC;            ///< temporal CABAC state storage for RD computation
  TEncBinCABAC            m_cRDGoOnBinCoderCABAC;         ///< going on bin coder CABAC for RD stage
#endif

  // quality control
  TEncPreanalyzer         m_cPreanalyzer;                 ///< image characteristics analyzer for TM5-step3-like adaptive QP

  TEncRateCtrl            m_cRateCtrl;                    ///< Rate control class

protected:
  Void  xGetNewPicBuffer  ( TComPic*& rpcPic, Int ppsId ); ///< get picture buffer which will be processed. If ppsId<0, then the ppsMap will be queried for the first match.
  Void  xInitVPS          (TComVPS &vps, const TComSPS &sps); ///< initialize VPS from encoder options
  Void  xInitSPS          (TComSPS &sps);                 ///< initialize SPS from encoder options
  Void  xInitPPS          (TComPPS &pps, const TComSPS &sps); ///< initialize PPS from encoder options
  Void  xInitScalingLists (TComSPS &sps, TComPPS &pps);   ///< initialize scaling lists
  Void  xInitHrdParameters(TComSPS &sps);                 ///< initialize HRD parameters

  Void  xInitPPSforTiles  (TComPPS &pps);
  Void  xInitRPS          (TComSPS &sps, Bool isFieldCoding);           ///< initialize PPS from encoder options

public:
  TEncTop();
  virtual ~TEncTop();

  Void      create          ();
  Void      destroy         ();
  Void      init            (Bool isFieldCoding);
  Void      deletePicBuffer ();

  Void      CRC_Seq_Init();

  // -------------------------------------------------------------------------------------------------------------------
  // member access functions
  // -------------------------------------------------------------------------------------------------------------------

  TComList<TComPic*>*     getListPic            () { return  &m_cListPic;             }
  TEncSearch*             getPredSearch         () { return  &m_cSearch;              }

  TComTrQuant*            getTrQuant            () { return  &m_cTrQuant;             }
  TComLoopFilter*         getLoopFilter         () { return  &m_cLoopFilter;          }
  TEncSampleAdaptiveOffset* getSAO              () { return  &m_cEncSAO;              }
  TEncGOP*                getGOPEncoder         () { return  &m_cGOPEncoder;          }
  TEncSlice*              getSliceEncoder       () { return  &m_cSliceEncoder;        }
  TEncCu*                 getCuEncoder          () { return  &m_cCuEncoder;           }
  TEncEntropy*            getEntropyCoder       () { return  &m_cEntropyCoder;        }
  TEncCavlc*              getCavlcCoder         () { return  &m_cCavlcCoder;          }
  TEncSbac*               getSbacCoder          () { return  &m_cSbacCoder;           }
  TEncBinCABAC*           getBinCABAC           () { return  &m_cBinCoderCABAC;       }

  TComRdCost*             getRdCost             () { return  &m_cRdCost;              }
  TEncSbac***             getRDSbacCoder        () { return  m_pppcRDSbacCoder;       }
  TEncSbac*               getRDGoOnSbacCoder    () { return  &m_cRDGoOnSbacCoder;     }
  TEncRateCtrl*           getRateCtrl           () { return &m_cRateCtrl;             }

  Int*					  getPicRcvd			() { return  &m_iNumPicRcvd; }	// Derek
  UInt*					  getAllPicCoded		() { return  &m_uiNumAllPicCoded; }

  Void selectReferencePictureSet(TComSlice* slice, Int POCCurr, Int GOPid );
  Int getReferencePictureSetIdxForSOP(Int POCCurr, Int GOPid );

#if JCTVC_Y0038_PARAMS
  Void                   setParamSetChanged(Int spsId, Int ppsId);
#endif
  Bool                   PPSNeedsWriting(Int ppsId);
  Bool                   SPSNeedsWriting(Int spsId);

  // -------------------------------------------------------------------------------------------------------------------
  // encoder function
  // -------------------------------------------------------------------------------------------------------------------

  /// encode several number of pictures until end-of-sequence
  Void encode( Bool bEos,
               TComPicYuv* pcPicYuvOrg,
               TComPicYuv* pcPicYuvTrueOrg, const InputColourSpaceConversion snrCSC, // used for SNR calculations. Picture in original colour space.
               TComList<TComPicYuv*>& rcListPicYuvRecOut,
               std::list<AccessUnit>& accessUnitsOut, Int& iNumEncoded );

  /// encode several number of pictures until end-of-sequence
  Void encode( Bool bEos, TComPicYuv* pcPicYuvOrg,
               TComPicYuv* pcPicYuvTrueOrg, const InputColourSpaceConversion snrCSC, // used for SNR calculations. Picture in original colour space.
               TComList<TComPicYuv*>& rcListPicYuvRecOut,
               std::list<AccessUnit>& accessUnitsOut, Int& iNumEncoded, Bool isTff);

  TEncAnalyze::OutputLogControl getOutputLogControl() const
  {
    TEncAnalyze::OutputLogControl outputLogCtrl;
    outputLogCtrl.printFrameMSE=m_printFrameMSE;
    outputLogCtrl.printMSEBasedSNR=m_printMSEBasedSequencePSNR;
#if JVET_F0064_MSSSIM
    outputLogCtrl.printMSSSIM=m_printMSSSIM;
#endif
    outputLogCtrl.printSequenceMSE=m_printSequenceMSE;
#if JCTVC_Y0037_XPSNR
    outputLogCtrl.printXPSNR=m_bXPSNREnableFlag;
#endif
    outputLogCtrl.printHexPerPOCPSNRs=m_printHexPsnr;
    return outputLogCtrl;
  }

  Void printSummary(Bool isField)
  {
    m_cGOPEncoder.printOutSummary (m_uiNumAllPicCoded, isField, getOutputLogControl(), m_spsMap.getFirstPS()->getBitDepths());
  }

};


#define CUSTOM_RC       1
#define QUICK_RC        0   // ADAPTIVE_QP_SELECTION=0 TransformSkipFast=1 FastDeltaQP=1 ECU=1 CFM=1
#define CTU_INFO_CSV    1
#define CTU_WEIGHT      1   // 1=high-pass, 2=luma/motion, 3=variance

typedef struct enc_frame_stat
{
    int frame_POC;
    int frame_SliceQP;
    int frame_AvgQP;
    int frame_maxQP;
    int frame_minQP;
    int slice_type;

    // pre-status
    //int estframeBits;
    //float m_nonSkipBits8x8;
    int m_avgGRD;
    int m_avgHAD;
    int m_avgSAD;
    int m_avgDHAD;
    //int m_avgDVAR;
    int m_avgDGR;
    UInt64 m_TotalSAD;

    // post-status
    int frameBits;   // without 0x02 0x03 and CABAC_ZERO_WORD, that is, only RBSP size
    int frameBitsHeader;
    int frameBitsData;
    int realSkip8x8Num;
    int m_TotalNonZeroNum_I;    // Only Y
    int m_TotalNonZeroNum_P;
    int m_TotalDistortion[3];
    //int m_MotionCnts[3];
    //int m_MotionBits[3];
    //int m_MotionSAD[3];
    //int m_MotionDGR[3];
    float m_encPSNR[4];
    float m_encSSIM[4];
    float m_encMSE[4];
    //int m_encMAD[4];

    int cpb_overflow;
    int cpb_state;
    int cpb_vbuf;
    int est_vbuf;
    int gop_left_bits;
    int target_bits;
    double pic_lambda;
    double pic_HAD_alpha;
    double pic_HAD_beta;

    double blk_weight_c[2];
}EncFrameStat;

typedef struct enc_ctu_text
{
    int m_TotalBits;
    int m_ctuHeaderBits;    // includes YUV, in raster scan
    int m_ctuDataBits;
    
    //int m_QP;
    int m_avgQP;
    int m_minQP;
    int m_maxQP;
    int m_CtuType[3]; // 0=skip, 1=inter, 2=intra
    float m_motion;   // 0=still, 1=small motion, 2=large motion
    int m_RdCost;
    int m_Dist[3];  // 0-Y, 1-UV, 2-YUV
    float m_PSNR[3];  // 0-Y, 1-U, 2-V
    int m_NonZeroCoeffNum;
    //int m_estSkip;

    int m_ctuGRD;  // intra-frame
    int m_ctuHAD;  // intra-frame
    int m_ctuSAD;  // inter-frame
    int m_ctuDHAD;  // inter-frame
    int m_ctuDVAR; // inter-frame
    int m_ctuDGR;  // inter-frame

    int m_ctuAvgLuma;
    int m_ctuAvgVAR;
    double motion_weight;
    double luma_weight;
    double lumo_weight; // V1
    double hipa_weight; // V2
    double var_weight;  // V3
    double m_CtuLambda;

    //int m_TargetBits;
    //int m_CPMode;
    int m_CtuLambdaIdx;
    int m_CtuLeftBits;
}EncCtuText;

typedef struct enc_top_stat
{
    vector<EncFrameStat>    FrameStat;
    int m_isTrial;
    int m_encFrameNum;
    TComPic* curPic;
    TComPic* refPic;

    EncCtuText *CtuText[2]; // cur and ref

    int m_maxQP;
    int m_minQP;
    int m_HeaderBits;  // temp variable
    int m_DataBits;    // temp variable
    int m_BitCnts;     // temp variable
    int m_CoeffNum_I;    // temp variable
    int m_CoeffNum_P;    // temp variable

    double HAD_ALPHA[52];
    double HAD_BETA[52];
}EncTopStat;
extern EncTopStat encTopStat;

typedef struct enc_top_crc
{
    // Seq
    int encFrameNum;
    int CRCWidth;
    int CRCHeight;
    int CRCPicPixels;
    int CRCCtuSize;

    int CRCEn;
    int CRCCtuEn;
    int CRCCpbEn;
    int CRCTargetBit;
    int CRCTargetAvgBit;
    int CRCInitQP;

    int CRCFrameRate;
    int CRCTotalFrames;
    int CRCFramesLeft;
    int CRCPrevPQP;
    double CRCHadScale;

    int CRCMaxCpbSize;
    int CRCCpbSize;
    int CRCCpbState;
    int CRCBufferingRate;

    // lambda table
    int CRCLambdaIdx;
    int CRCLambdaIdxAccP;
    int CRCLambdaIdxAccN;
    double CRCLambdaTBL[52*4];
    double CRCPicConstI;
    double CRCPicConstP;
    double CRCPicConstPAvg;
    double CRCPicConstPSum;
    double CRCPicConstPSlot[32];    // CONSTP_WIN
    int CRCPicConstPBitsAvg;
    int CRCPicConstPBitsSum;
    int CRCPicConstPBitsSlot[32];   // CONSTP_WIN
    int CRCPicConstPPrvIdx;
    int CRCPicConstPCnt;

    // Gop
    int CRCGopSize;
    int CRCGopLeftBits;
    int CRCGopDeltaBits;
    int CRCGopTargetBits;
    int CRCGopTargetVBuf;
    int CRCGopAvgVBuf;

    // Pic
    int CRCPicIdealBits;
    int CRCPicIdealBits2;
    int CRCPicTargetBits;

    int CRCPicAvgIBits;
    int CRCPicAvgIBitsSum;
    int CRCPicAvgIBitsSlot[4];   // SMOOTH_I_WIN
    int CRCPicAvgIPrvIdx;
    int CRCPicAvgICnt;

    int CRCPicAvgLeftBits;
    int CRCPicAvgLeftBitsSum;
    int CRCPicAvgLeftBitsSlot[4];   // LEFTBITS_WIN
    int CRCPicAvgLeftBitsPrvIdx;
    int CRCPicAvgLeftBitsCnt;

    int CRCPicQP;
    int CRCPicCompMode;
    double CRCPicLambda;

    // HAD
    int CRCPicEstHADIBits;
    int CRCPicEstHADIQP;

    // SAD
    int CRCPicSAD;
    int CRCPicPrvSAD;
    //int CRCPicPrvLambdaIdx;
    int CRCPicPrvLambdaOffset;
    int CRCPicSADLambdaShift;

    int CRCPicAvgSAD;
    int CRCPicAvgSADSum;
    int CRCPicAvgSADSlot[16];   // SAD_WIN

    int CRCPicAvgLambdaIdx;
    int CRCPicAvgLambdaIdxSum;
    int CRCPicAvgLambdaIdxSlot[16];   // SAD_WIN

    double CRCPicAvgSADAlpha;
    double CRCPicAvgSADAlphaSum;
    double CRCPicAvgSADAlphaSlot[16];   // SAD_WIN

    double CRCPicAvgSADBeta;
    double CRCPicAvgSADBetaSum;
    double CRCPicAvgSADBetaSlot[16];   // SAD_WIN

    int CRCPicAvgSADConst;
    int CRCPicAvgSADConstSum;
    int CRCPicAvgSADConstSlot[16];   // SAD_WIN

    double CRCPicSADAlpha;
    double CRCPicSADBeta;
    int CRCPicSADConst;

    int CRCPicAvgSADPrvIdx;
    int CRCPicAvgSADCnt;

    // CTU
    //int CRCCtuTargetBits;
    //int CRCCtuActualBits;
    int CRCCtuAccumBits;
    //int CRCCtuCompMode;
    int CRCCtuLambdaIdx;
    int CRCCtuQP;
    int CRCQPC_delta[2];
    double CRCCtuLambda;
    double CRCCtuLambdaScale;
    int CRCCtuAvgLambdaIdx;
    //int CRCCtuAvgQP;
    int CRCCtuValidCnt;
    int CRCCtuLambdaIdxOffset;

    // VBR
    int CRCPicVBRMinBits;
    int CRCPicVBRMaxBits;
    int CRCPicVBRMinQP;
    int CRCPicVBRMaxQP;    
    int CRCPicPrvIBits;
    int CRCPicVBRHits;
    int CRCPicEstActBits;

    //int CRCPicVBRAvgPBits;
    //int CRCPicVBRPBitsSum;
    //int CRCPicVBRPBitsSlot[16];   // VBR_P_WIN
    //int CRCPicVBRPrvPIdx;
    //int CRCPicVBRPCnt;
    
}EncTopCRC;
extern EncTopCRC encTopCRC;

#if CHK_RDO_HAD
extern int g_RDOTotalNum;
extern int g_HADEqRDONum;
extern double g_RatioSum;
extern FILE *fp_RDOHAD;
#endif

void Enc_Init_Text(int width, int height, int ctu_size);
void Enc_DeInit_Text();
void Analyze_Pictures(TComPic* currPic, int intra_period);
void Estimate_Intra_Bits();
int Estimate_Intra_QP(int target_bits);
void update_parameters(int actual_bits);
void Print_CTU_Text(TComPic* pcPic, int intra_period);
void Print_Frame_Stat();

int CRC_Init_GOP();
void CRC_Init_PIC(int isISlice);
void CRC_Update_Params(int isISlice);

int find_lambda_idx(int Idx, const double estLambda_cur);
void CRC_Estimate_CTU_Lambda(TComPic* pcPic, int totalCtuInSlice, int CtuLeftNum);
//! \}

#endif // __TENCTOP__

