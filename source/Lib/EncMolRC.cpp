
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include "../App/TAppEncoder/TAppEncTop.h"
#include "TLibEncoder/TEncTop.h"

#if ORG_TBL
const double HAD_ALPHA_POW[52] =
{
    564927, 522095, 473260, 424316, 337334, 280559,
    223524, 171344, 120886, 87352,  84368,  40963,
    41124,  20339,  19092,  9448.8, 8697.2, 4741.7,
    4292.1, 2430.3, 2248.5, 1300.5, 1280.4, 782.23,
    805.06, 507,    487.4,  321.72, 296.15, 206.38,
    182.07, 130.2,  100.26, 72.597, 52.002, 41.664,
    28.593, 25.991, 16.891, 18.498, 13.208, 16.591,
    12.616, 17.018, 14.001, 17.862, 26.361, 28.714,
    32.467, 37.261, 46.32, 57.591
};
const double HAD_BETA_POW[52] =
{
    0.4295, 0.4363, 0.446,  0.4573, 0.4827, 0.5016,
    0.5269, 0.5559, 0.5954, 0.6326, 0.6269, 0.718,
    0.7101, 0.7982, 0.7983, 0.8864, 0.8889, 0.9631,
    0.9677, 1.0366, 1.0359, 1.1019, 1.0915, 1.1498,
    1.1326, 1.186,  1.1758, 1.2223, 1.2178, 1.2554,
    1.2578, 1.291,  1.3089, 1.3398, 1.3673, 1.3833,
    1.4174, 1.4145, 1.4538, 1.4247, 1.4517, 1.4023,
    1.4211, 1.3612, 1.366,  1.3122, 1.233,  1.1985,
    1.1595, 1.1162, 1.0598, 1.0046
};

const double HAD_ALPHA_LIN[52] =
{
    3901.3, 3821.9, 3757.6, 3722,   3708.3, 3649,
    3615.5, 3573.5, 3522.9, 3507.2, 3222.1, 3193.1,
    3124.3, 3051.4, 2969.1, 2882.1, 2804.6, 2673.7,
    2603.7, 2462.2, 2362.6, 2229.8, 2126.1, 1995.5,
    1899.2, 1778.4, 1683,   1575.9, 1495.5, 1386,
    1326.5, 1222.2, 1147.7, 1054.2, 973.85, 883.71,
    813.16, 728.07, 658.55, 582.88, 521.41, 452.75,
    401.75, 344.19, 296.02, 249.59, 202.91, 168.11,
    139.16, 112.34, 87.95,  69.03
};
const double HAD_BETA_LIN[52] =
{
    7000000, 7000000, 6000000, 6000000, 6000000, 5000000,
    5000000, 4000000, 4000000, 3000000, 3000000, 3000000,
    2000000, 2000000, 2000000, 1000000, 1000000, 988390,
    833820,  663841,  531622,  405481,  303799,  219410,
    137655,  73007,   -15100,  -60494,  -136828, -158859,
    -221256, -226871, -276999, -270980, -293095, -270836,
    -279282, -248180, -247854, -211995, -202472, -166341,
    -154809, -123277, -109403, -85019,  -63142,  -48807,
    -36680,  -25767,  -15573,  -8012.2
};
#endif

#if TBL2
double HAD_ALPHA_LIN[52] =
{
    3901.3, 3821.9, 3757.6, 3722.0, 3708.3, 3649.0,
    3615.5, 3573.5, 3522.9, 3507.2, 3222.1, 3193.1,
    3124.3, 3051.4, 2969.1, 3826.5, 2804.6, 2673.7,
    2603.7, 2462.2, 2362.6, 2229.8, 2126.1, 1995.5,
    1899.2, 1778.4, 1683.0, 1575.9, 1495.5, 1350.2,
    1226.0, 1105.5, 976.00, 875.80, 768.59, 688.01,
    605.93, 539.40, 467.08, 416.45, 361.91, 319.11,
    277.19, 243.71, 206.18, 178.31, 149.22, 125.62,
    105.86, 87.770, 71.091, 57.218
};

double HAD_BETA_LIN[52] =
{
    7000000, 7000000, 6000000, 6000000, 6000000, 5000000,
    5000000, 4000000, 4000000, 3000000, 3000000, 3000000,
    2000000, 2000000, 2000000, 1000000, 1000000, 988390,
    833820,  663841,  531622,  405481,  303799,  219410,
    137655,  73007,  -15100,  -60494,  -136828, -120000,
   -110000, -100000, -87000,  -77000,  -66000,  -58000,
   -50000,  -43000,  -36000,  -31000,  -26000,  -21000,
   -17000,  -14000,  -10000,  -7500,   -4000,   -2000,
    0, 1300, 3000, 5000
};

int MIN_HAD_LIN[52] =
{
    1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1,
    1, 1, 50, 100, 100, 100,
    100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100,
    100, 100, 100, 100
};

int MAX_HAD_LIN[52] =
{
    1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1,
    1, 1, 1, 70, 80, 100,
    150, 250, 300, 400, 450, 500,
    600, 650, 750, 800, 850, 950,
    1050, 1150, 1300, 1450, 1600, 1800,
    2000, 2000, 2500, 3000, 3500, 4000,
    4000, 5000, 6000, 7000, 8000, 10000,
    10000, 10000, 10000, 10000
};
#endif

double HAD_ALPHA_LIN[52] =
{
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 6599.1, 6069.2,
    5725.8, 5245.8, 4893.7, 4474.0, 4167.3, 3802.7,
    3535.6, 3224.8, 2955.0, 2694.9, 2459.3, 2242.1,
    2040.2, 1845.5, 1669.4, 1520.2, 1371.8, 1239.9,
    1226.0, 1105.5, 976.00, 875.80, 768.59, 688.01,
    605.93, 539.40, 467.08, 416.45, 361.91, 319.11,
    277.19, 243.71, 206.18, 178.31, 149.22, 125.62,
    105.86, 87.770, 71.091, 57.218
};

double HAD_BETA_LIN[52] =
{
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    -110000, -100000, -87000, -77000, -66000, -58000,
    -50000,  -43000,  -36000, -31000, -26000, -21000,
    -17000,  -14000,  -10000, -7500,  -4000,  -2000,
    0, 1300, 3000, 5000
};

int MIN_HAD_LIN[52] =
{
    1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1,
    1, 5, 5, 5, 5, 5,
    100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100,
    100, 100, 100, 100
};

int MAX_HAD_LIN[52] =
{
    1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 30, 40,
    50, 60, 70, 80, 90, 100,
    120, 140, 160, 200, 250, 300,
    400, 550, 650, 750, 850, 1000,
    1200, 1350, 1500, 1700, 2000, 2200,
    2500, 2800, 3200, 3500, 4000, 4500,
    5200, 6000, 7000, 8000, 10000, 10000,
    10000, 10000, 10000, 10000
};

#define RESHAPE_BITS    0

#define FRM_NUM         50
#define DEFAULT_RATE    4000000
#define MIN_BITS_I      400000  // bit
#define MIN_BITS_P      20000   // bit
#define MAX_BITS_I      (DEFAULT_RATE - (FRM_NUM-1)*MIN_BITS_P)
#define MAX_CRC_QP      51
#define MIN_CRC_QP      10
#define MAX_CRC_DQP_P   2
#define MAX_CRC_DQP_I_P   2
#define MAX_CRC_DQP_I_N   2
#define SMOOTH_I_WIN    4
#define SMOOTH_I_WT     0.7
#define LEFTBITS_WIN    4
#define LEFTBITS_WT     0.7
#define CONSTP_WIN      32

#define CONVERGE_FACTOR 0.25
#define COMPRESS_FACTOR 0.1
#define I_TO_P_EST      1.25    // smaller, less I to P drop
#define P_TO_I_EST      1.0    // larger, larger P to I raise
#define LAMBDA_NUMQP    4

#define VBR_EN          0   // set CBR_SAD_EN=0
#define VBR_QUALITY     5   // 1 to 10, higher is better
#define SAD_WIN         8

#define CBR_SAD_EN      (0 && !VBR_EN)
#define VBR_SAD_EN      VBR_EN
#define PIC_RATE        (DEFAULT_RATE / FRM_NUM)
#define VBR_MAX_QP      40
#define VBR_MIN_QP      15
#define VBR_MAX_RATE    int(2*PIC_RATE)
#define VBR_MIN_RATE    int(0.1*PIC_RATE)
//#define VBR_P_WIN       16

#define SAD_EN          (CBR_SAD_EN || VBR_SAD_EN)

Void TAppEncTop::CRC_Seq_Init()
{
    encTopCRC.CRCEn = m_CRCEnable;
    encTopCRC.CRCCtuEn = m_CRCCtuEnable;
    encTopCRC.CRCCpbEn = m_CRCCpbEnable;
    encTopCRC.CRCTargetBit = m_CRCTargetBitrate;
    encTopCRC.CRCInitQP = encTopCRC.CRCPrevPQP = m_CRCInitialQP;
    encTopCRC.CRCFrameRate = m_iFrameRate;
    encTopCRC.CRCTotalFrames = encTopCRC.CRCFramesLeft = m_framesToBeEncoded;
    encTopCRC.CRCGopSize = m_iIntraPeriod;

    const UInt uiLevelIdx = (m_level / 10) + (UInt)((m_level % 10) / 3);    // (m_level / 30)*3 + ((m_level % 10) / 3);
    const UInt g_uiMaxCpbSize[2][21] =
    {
        //         LEVEL1,        LEVEL2,LEVEL2_1,     LEVEL3, LEVEL3_1,      LEVEL4, LEVEL4_1,       LEVEL5,  LEVEL5_1,  LEVEL5_2,    LEVEL6,  LEVEL6_1,  LEVEL6_2
        { 0, 0, 0, 350000, 0, 0, 1500000, 3000000, 0, 6000000, 10000000, 0, 12000000, 20000000, 0,  25000000,  40000000,  60000000,  60000000, 120000000, 240000000 },
        { 0, 0, 0,      0, 0, 0,       0,       0, 0,       0,        0, 0, 30000000, 50000000, 0, 100000000, 160000000, 240000000, 240000000, 480000000, 800000000 }
    };
    encTopCRC.CRCMaxCpbSize = g_uiMaxCpbSize[m_levelTier][uiLevelIdx];
    encTopCRC.CRCWidth = m_iSourceWidth;
    encTopCRC.CRCHeight = m_iSourceHeight;
    encTopCRC.CRCPicPixels = encTopCRC.CRCWidth * encTopCRC.CRCHeight;
    encTopCRC.CRCCtuSize = m_uiMaxCUWidth;
    const Int defaultCPB = encTopCRC.CRCTargetBit; // 1 sec // Int(encTopCRC.CRCWidth * encTopCRC.CRCHeight * 0.8 + 0.5);
    encTopCRC.CRCCpbSize = min(defaultCPB, encTopCRC.CRCMaxCpbSize);
    encTopCRC.CRCCpbState = Int(encTopCRC.CRCCpbSize * 0.9 + 0.5);
    encTopCRC.CRCBufferingRate = (encTopCRC.CRCTargetBit + encTopCRC.CRCFrameRate - 1) / encTopCRC.CRCFrameRate;
    encTopCRC.CRCGopTargetBits = 0;
    encTopCRC.CRCGopDeltaBits = 0;
    encTopCRC.CRCGopLeftBits = 0;
    encTopCRC.CRCPicAvgIBits = 0;
    encTopCRC.CRCPicCompMode = 0;
    encTopCRC.CRCCtuLambdaScale = 1;
    encTopCRC.CRCCtuAvgLambdaIdx = 0;
    //encTopCRC.CRCCtuAvgQP = 0;
    encTopCRC.CRCCtuValidCnt = 0;
    encTopCRC.CRCCtuLambdaIdxOffset = 0;

    // reset Sliding Window for ConstP
    encTopCRC.CRCPicConstPCnt = 0;
    encTopCRC.CRCPicConstPPrvIdx = 0;
    encTopCRC.CRCPicConstPSum = 0;
    encTopCRC.CRCPicConstPBitsSum = 0;

#if SAD_EN
    // reset Sliding Window for Average SAD
    encTopCRC.CRCPicAvgSAD = 0;
    encTopCRC.CRCPicAvgSADCnt = 0;
    encTopCRC.CRCPicAvgSADPrvIdx = 0;
    encTopCRC.CRCPicAvgSADSum = 0;
    encTopCRC.CRCPicAvgLambdaIdxSum = 0;
    encTopCRC.CRCPicAvgSADAlphaSum = 0;
    encTopCRC.CRCPicAvgSADBetaSum = 0;
    encTopCRC.CRCPicAvgSADConstSum = 0;
    //encTopCRC.CRCPicSADBits = 0;
    encTopCRC.CRCPicSADAlpha = 1000;
    encTopCRC.CRCPicSADBeta = 1000;
    encTopCRC.CRCPicSADConst = 0;
    encTopCRC.CRCPicSADLambdaShift = 0;
    encTopCRC.CRCPicPrvLambdaOffset = 0;
#endif

    // rescale HAD coeffs ?
    encTopCRC.CRCHadScale = (double)encTopCRC.CRCWidth * encTopCRC.CRCHeight / 1920 / 1088;
    if (encTopCRC.CRCWidth != 1920 && encTopCRC.CRCHeight != 1088)
    {
        for (int i = 0; i < 52; i++)
        {
            HAD_ALPHA_LIN[i] *= encTopCRC.CRCHadScale;
            HAD_BETA_LIN[i] *= encTopCRC.CRCHadScale;
        }
    }

#if VBR_EN
    encTopCRC.CRCPicVBRMinQP = VBR_MIN_QP;
    encTopCRC.CRCPicVBRMaxQP = VBR_MAX_QP;
    encTopCRC.CRCPicVBRMinBits = int(VBR_MIN_RATE * encTopCRC.CRCHadScale);
    encTopCRC.CRCPicVBRMaxBits = int(VBR_MAX_RATE * encTopCRC.CRCHadScale);
    encTopCRC.CRCPicVBRHits = 1;

    // reset Sliding Window for VBR PBits
    //encTopCRC.CRCPicVBRPCnt = 0;
    //encTopCRC.CRCPicVBRPrvPIdx = 0;
    //encTopCRC.CRCPicVBRPBitsSum = 0;
#endif

    // initial lambda table
    for (int QP = 0; QP <= 51; QP++)
    {
        Double lambda = exp((((Double)QP - 13.7122) / 4.2005));
        Double lambda_diff = 1.0613235; // for LAMBDA_NUMQP, lambda^2 for QPs is 1.6098383876173012
        for (int idx = 0; idx < LAMBDA_NUMQP; idx++)
        {
            encTopCRC.CRCLambdaTBL[QP * LAMBDA_NUMQP + idx] = lambda * lambda_diff;
            lambda = encTopCRC.CRCLambdaTBL[QP * LAMBDA_NUMQP + idx];
        }
    }

    // test
    /*double a[13];
    for (int idx = 1; idx < 13; idx++)
    {
        a[idx] = pow(1.12640, idx);
    }
    // test
    double lambda_scale0[52];
    for (int idx = 1; idx < 52; idx++)
    {
        Double lambda_c = encTopCRC.CRCLambdaTBL[idx*LAMBDA_NUMQP];
        Double lambda_c2 = lambda_c * lambda_c;
        Double lambda_p = encTopCRC.CRCLambdaTBL[idx*LAMBDA_NUMQP - LAMBDA_NUMQP];
        Double lambda_p2 = lambda_p * lambda_p;

        lambda_scale0[idx] = lambda_c2 / lambda_p2;
    }
    // test, 1.1264
    double lambda_scale[52 * LAMBDA_NUMQP];
    for (int idx = 1; idx < 52 * LAMBDA_NUMQP; idx++)
    {
        Double lambda_c = encTopCRC.CRCLambdaTBL[idx];
        Double lambda_c2 = lambda_c * lambda_c;
        Double lambda_p = encTopCRC.CRCLambdaTBL[idx-1];
        Double lambda_p2 = lambda_p * lambda_p;

        lambda_scale[idx] = lambda_c2 / lambda_p2;
    }
    lambda_scale[0] = 0;*/
}

int CRC_Init_GOP()
{
    encTopCRC.CRCTargetAvgBit = encTopCRC.CRCTargetBit / encTopCRC.CRCFrameRate; // reset target bits
    encTopCRC.CRCGopAvgVBuf = encTopCRC.CRCGopDeltaBits / encTopCRC.CRCGopSize; //g_RCSmoothWindowSize, Derek, can be other time period ?
    encTopCRC.CRCPicIdealBits = encTopCRC.CRCTargetAvgBit - encTopCRC.CRCGopAvgVBuf;
    int targetBits = encTopCRC.CRCPicIdealBits * encTopCRC.CRCGopSize;

    encTopCRC.CRCGopTargetVBuf = encTopCRC.CRCGopDeltaBits;

    if (targetBits < int(MIN_BITS_P * FRM_NUM * encTopCRC.CRCHadScale))
    {
        targetBits = int(MIN_BITS_P * FRM_NUM * encTopCRC.CRCHadScale);   // at least allocate X bits for one GOP
    }
    return targetBits;
}

int CRC_Estimate_IBits2(int est_QP, const int estIBits_cur)
{
    const int frameno = encTopStat.m_encFrameNum;
    const int HAD = encTopStat.FrameStat[frameno].m_avgHAD;

    double alpha = encTopStat.HAD_ALPHA[est_QP];
    double beta = encTopStat.HAD_BETA[est_QP];
    int estIBits_opt = int(alpha * HAD + beta);
    const int max_QP = VBR_EN ? encTopCRC.CRCPicVBRMaxQP : MAX_CRC_QP;
    const int min_QP = VBR_EN ? encTopCRC.CRCPicVBRMinQP : MIN_CRC_QP;

    // Find Optimal QP according to bits
    int delta_more, delta_less;
    if (estIBits_cur > estIBits_opt)
    {
        do
        {
            delta_more = estIBits_cur - estIBits_opt;
            est_QP--;
            if (est_QP < min_QP)
            {
                est_QP = min_QP;
                break;
            }

            alpha = encTopStat.HAD_ALPHA[est_QP];
            beta = encTopStat.HAD_BETA[est_QP];
            estIBits_opt = int(alpha * HAD + beta);
        } while (estIBits_cur > estIBits_opt);

        delta_less = estIBits_opt - estIBits_cur;
        if (delta_more < delta_less)
            est_QP++;
    }
    else    // (estIBits_cur < estIBits_opt)
    {
        do
        {
            delta_less = estIBits_opt - estIBits_cur;
            est_QP++;
            if (est_QP > max_QP)
            {
                est_QP = max_QP;
                break;
            }

            alpha = encTopStat.HAD_ALPHA[est_QP];
            beta = encTopStat.HAD_BETA[est_QP];
            estIBits_opt = int(alpha * HAD + beta);

        } while (estIBits_cur < estIBits_opt);

        delta_more = estIBits_cur - estIBits_opt;
        if (delta_less < delta_more)
            est_QP--;
    }
    return est_QP;
}

int CRC_Estimate_IBits(int est_QP)
{
    int QP = est_QP;
    const int frameno = encTopStat.m_encFrameNum;
    const int HAD = encTopStat.FrameStat[frameno].m_avgHAD;
#if VBR_EN
    const int max_Ibits = int(MAX_BITS_I * encTopCRC.CRCHadScale);
#else
    const int max_Ibits = encTopCRC.CRCGopTargetBits - int((encTopCRC.CRCGopSize - 1) * MIN_BITS_P * encTopCRC.CRCHadScale);
#endif
    const int min_Ibits = int(MIN_BITS_I * encTopCRC.CRCHadScale);

    if (frameno == 0)
    {
        // Find Optimal QP First according to HAD
        int MAX_HAD = MAX_HAD_LIN[est_QP];
        if (HAD <= MAX_HAD)
        {
            do
            {
                est_QP--;
                MAX_HAD = MAX_HAD_LIN[est_QP];
            } while (HAD < MAX_HAD);
            est_QP++;
        }
        else
        {
            do
            {
                est_QP++;
                MAX_HAD = MAX_HAD_LIN[est_QP];
            } while (HAD > MAX_HAD);
        }

        QP = est_QP;
    }
    else
    {
        // bits from previous GOPs
        const int bitsFromPrvGOP = encTopCRC.CRCPicAvgIBits;
        // bits from previous P
        const int bitsFromP = int(encTopCRC.CRCPicPixels * encTopCRC.CRCPicConstI / encTopCRC.CRCPicLambda / encTopCRC.CRCPicLambda);
        const int estBits = int(SMOOTH_I_WT * bitsFromPrvGOP + (1 - SMOOTH_I_WT) * bitsFromP);
        int old_qp = est_QP;
        est_QP = CRC_Estimate_IBits2(est_QP, estBits);
        printf("RC: I Slice (BitsP, BitsGOP, QP, oldQP) = (%d, %d, %d, %d)\n", bitsFromP, bitsFromPrvGOP, est_QP, old_qp);
        {
            const int bits0 = int(encTopStat.HAD_ALPHA[est_QP - 1] * HAD + encTopStat.HAD_BETA[QP - 1]);
            const int bits1 = int(encTopStat.HAD_ALPHA[est_QP + 1] * HAD + encTopStat.HAD_BETA[QP + 1]);
            printf("RC: I Slice (estBits, QP+1, QP-1) = (%d, %d, %d)\n", estBits, bits1, bits0);
        }
        QP = Clip3(QP - MAX_CRC_DQP_I_N, QP + MAX_CRC_DQP_I_P, est_QP);   // prevent IP flickering, P to I uses smaller ?
    }

    QP = Clip3(MIN_CRC_QP, MAX_CRC_QP, QP);
#if VBR_EN
    QP = Clip3(encTopCRC.CRCPicVBRMinQP, encTopCRC.CRCPicVBRMaxQP, QP);
#endif
    // estimate bits
    const double alpha = encTopStat.HAD_ALPHA[QP];
    const double beta = encTopStat.HAD_BETA[QP];
    const int bits = int(alpha * HAD + beta);
    encTopCRC.CRCPicEstHADIQP = QP;
    return Clip3(min_Ibits, max_Ibits, bits);
}

int CRC_VBR_QPADJUST(int inbits)
{
    const int frameno = encTopStat.m_encFrameNum;
    int outbits = inbits;
    // Calculate SAD lambda offsets
    if (encTopCRC.CRCPicAvgSADCnt != 0)
    {
        encTopCRC.CRCPicAvgSAD = encTopCRC.CRCPicAvgSADSum / encTopCRC.CRCPicAvgSADCnt;
        encTopCRC.CRCPicAvgLambdaIdx = encTopCRC.CRCPicAvgLambdaIdxSum / encTopCRC.CRCPicAvgSADCnt;
        encTopCRC.CRCPicAvgSADAlpha = encTopCRC.CRCPicAvgSADAlphaSum / encTopCRC.CRCPicAvgSADCnt;
        encTopCRC.CRCPicAvgSADBeta = encTopCRC.CRCPicAvgSADBetaSum / encTopCRC.CRCPicAvgSADCnt;
        encTopCRC.CRCPicAvgSADConst = encTopCRC.CRCPicAvgSADConstSum / encTopCRC.CRCPicAvgSADCnt;

        int AvgSADBitsA = int(encTopCRC.CRCPicAvgSADAlpha * encTopCRC.CRCPicAvgSAD / encTopCRC.CRCPicAvgLambdaIdx + encTopCRC.CRCPicAvgSADBeta);
        int AvgPrvCurSAD = int(0.75 * encTopCRC.CRCPicSAD + 0.25 * encTopCRC.CRCPicPrvSAD);
        double estLambdaIdxAvgA = (encTopCRC.CRCPicAvgSADAlpha * AvgPrvCurSAD) / (AvgSADBitsA - encTopCRC.CRCPicAvgSADBeta);
        int lambdaOffsetAvgA = (int(estLambdaIdxAvgA + 0.5) - encTopCRC.CRCLambdaIdx);
        int sign = (lambdaOffsetAvgA >= 0) ? 1 : -1;
        int val = abs(lambdaOffsetAvgA);

        int diff = 0;
        if (val > 16)
        {
            if (sign > 0)    // large SAD
            {
                diff = ((val >> 4) + (val >> 5) + 1) >> 1;
                outbits = int(inbits * (1 + 0.05 * min(16, diff)));
                encTopCRC.CRCPicCompMode = 2;   // raise bits
            }
            else if (sign < 0)    // small SAD
            {
                diff = (val >> 4);
                //diff = ((val >> 3) + (val >> 4) + 1) >> 1;
                outbits = int(inbits * (1 - 0.05 * min(16, diff)));
                encTopCRC.CRCPicCompMode = 1;   // lower bits
            }
            encTopCRC.CRCPicVBRHits = 1;
        }
        else
        {
            int amount = min(encTopCRC.CRCPicVBRHits, 10-VBR_QUALITY) * 8;
            outbits = int(inbits * (100 - amount) / 100);
            encTopCRC.CRCPicVBRHits++;
        }

        // Est Prv Actual Bits
        const int PrvActualBits = encTopStat.FrameStat[frameno - 1].frameBits;
        int CurSADBits = int(encTopCRC.CRCPicAvgSADAlpha * AvgPrvCurSAD / encTopCRC.CRCLambdaIdx + encTopCRC.CRCPicAvgSADBeta);
        int PrvSADBits = int(encTopCRC.CRCPicAvgSADAlpha * encTopCRC.CRCPicPrvSAD / encTopCRC.CRCLambdaIdx + encTopCRC.CRCPicAvgSADBeta);
        encTopCRC.CRCPicEstActBits = int((double)CurSADBits / PrvSADBits * PrvActualBits);

        printf("RC: SAD(CP=%d)(Cur, PrvC, AvgA)=(%d, %d, %d), (AvgBits, SADEstBits, diffVal, EstPrvBits)=(%d, %d, %d, %d)\n", encTopCRC.CRCPicCompMode,
               encTopCRC.CRCPicSAD, AvgPrvCurSAD, encTopCRC.CRCPicAvgSAD, inbits, outbits, diff, encTopCRC.CRCPicEstActBits);
    }

    int limit_SAD_min, limit_SAD_max, limit_SAD;
    limit_SAD = encTopCRC.CRCPicSAD;
    if (encTopCRC.CRCPicAvgSAD != 0)
    {
        limit_SAD_min = int(0.5 * encTopCRC.CRCPicAvgSAD);
        limit_SAD_max = int(1.5 * encTopCRC.CRCPicAvgSAD);
        limit_SAD = Clip3(limit_SAD_min, limit_SAD_max, limit_SAD);
    }

    // Sliding Window for Average SAD
    encTopCRC.CRCPicAvgSADSum += limit_SAD;
    encTopCRC.CRCPicAvgLambdaIdxSum += encTopCRC.CRCLambdaIdx;
    encTopCRC.CRCPicAvgSADAlphaSum += encTopCRC.CRCPicSADAlpha;
    encTopCRC.CRCPicAvgSADBetaSum += encTopCRC.CRCPicSADBeta;
    encTopCRC.CRCPicAvgSADConstSum += encTopCRC.CRCPicSADConst;
    encTopCRC.CRCPicAvgSADCnt++;
    if (encTopCRC.CRCPicAvgSADCnt > SAD_WIN)
    {
        encTopCRC.CRCPicAvgSADSum -= encTopCRC.CRCPicAvgSADSlot[encTopCRC.CRCPicAvgSADPrvIdx];
        encTopCRC.CRCPicAvgLambdaIdxSum -= encTopCRC.CRCPicAvgLambdaIdxSlot[encTopCRC.CRCPicAvgSADPrvIdx];
        encTopCRC.CRCPicAvgSADAlphaSum -= encTopCRC.CRCPicAvgSADAlphaSlot[encTopCRC.CRCPicAvgSADPrvIdx];
        encTopCRC.CRCPicAvgSADBetaSum -= encTopCRC.CRCPicAvgSADBetaSlot[encTopCRC.CRCPicAvgSADPrvIdx];
        encTopCRC.CRCPicAvgSADConstSum -= encTopCRC.CRCPicAvgSADConstSlot[encTopCRC.CRCPicAvgSADPrvIdx];
        encTopCRC.CRCPicAvgSADCnt--;
    }
    encTopCRC.CRCPicAvgSADSlot[encTopCRC.CRCPicAvgSADPrvIdx] = limit_SAD;
    encTopCRC.CRCPicAvgLambdaIdxSlot[encTopCRC.CRCPicAvgSADPrvIdx] = encTopCRC.CRCLambdaIdx;
    encTopCRC.CRCPicAvgSADAlphaSlot[encTopCRC.CRCPicAvgSADPrvIdx] = encTopCRC.CRCPicSADAlpha;
    encTopCRC.CRCPicAvgSADBetaSlot[encTopCRC.CRCPicAvgSADPrvIdx] = encTopCRC.CRCPicSADBeta;
    encTopCRC.CRCPicAvgSADConstSlot[encTopCRC.CRCPicAvgSADPrvIdx] = encTopCRC.CRCPicSADConst;
    encTopCRC.CRCPicAvgSADPrvIdx = (encTopCRC.CRCPicAvgSADPrvIdx + 1) % SAD_WIN;

    return outbits;
}

int CRC_Estimate_PBits()
{
    int bits;
    const int GopSize = encTopCRC.CRCGopSize;
#if VBR_EN
    const int GopRate = int(encTopCRC.CRCBufferingRate * GopSize * 0.8);
    const int GopTotalBits = encTopCRC.CRCPicPrvIBits + (GopSize - 1) * encTopCRC.CRCPicConstPBitsAvg;
    int avgTargetBits = (UInt64)encTopCRC.CRCPicConstPBitsAvg * GopRate / GopTotalBits;
        
    int estTargetBits = avgTargetBits;
    const int frameno = encTopStat.m_encFrameNum;
    encTopCRC.CRCPicSAD = encTopStat.FrameStat[frameno].m_avgSAD;
    if (encTopCRC.CRCPicSAD != 0)
    {
        estTargetBits = CRC_VBR_QPADJUST(avgTargetBits);
    }
    bits = Clip3(encTopCRC.CRCPicVBRMinBits, encTopCRC.CRCPicVBRMaxBits, estTargetBits);
    encTopCRC.CRCPicCompMode = 0;
#else
    const int frameno = encTopStat.m_encFrameNum;
    const int frm_pos = frameno % GopSize;
    const int left_frm = GopSize - frm_pos;

    const int ideal_bits = encTopCRC.CRCPicIdealBits2;
    int gop_left_bits = encTopCRC.CRCGopLeftBits / left_frm;  // Derek, use SAD to distribute bits?
    const int VBuf_diff_bits = encTopCRC.CRCGopTargetVBuf - encTopCRC.CRCGopDeltaBits;

    // use diff weight for pic position
    const double pos_weight = 0.5;// +(0.3 * frm_pos) / GopSize;

    // detect if compress more
    int target_bits = ideal_bits;
    int old_gop_left = gop_left_bits;
    if (gop_left_bits < int((1 - COMPRESS_FACTOR) * target_bits))
    {
        encTopCRC.CRCPicCompMode = 1;
        // clip gop_left_bits if unreasonable
        if (gop_left_bits < int(0.8 * target_bits))
            gop_left_bits = int(0.8 * target_bits);
    }
    else if (gop_left_bits > int((1 + COMPRESS_FACTOR) * target_bits))
    {
        encTopCRC.CRCPicCompMode = 2;
        // clip gop_left_bits if unreasonable
        if (gop_left_bits > int(1.2 * target_bits))
            gop_left_bits = int(1.2 * target_bits);
    }
    else
        encTopCRC.CRCPicCompMode = 0;

    // Sliding Window for Average LeftBits, appied when left_frm is too low
    int real_left_bits = gop_left_bits;
    if (left_frm <= LEFTBITS_WIN)
    {
        encTopCRC.CRCPicAvgLeftBits = encTopCRC.CRCPicAvgLeftBitsSum / encTopCRC.CRCPicAvgLeftBitsCnt;
        gop_left_bits = int(LEFTBITS_WT * encTopCRC.CRCPicAvgLeftBits + (1 - LEFTBITS_WT) * gop_left_bits);
    }
    if (left_frm <= 2 * LEFTBITS_WIN)
    {
        encTopCRC.CRCPicAvgLeftBitsSum += real_left_bits;
        encTopCRC.CRCPicAvgLeftBitsCnt++;
        if (encTopCRC.CRCPicAvgLeftBitsCnt > LEFTBITS_WIN)
        {
            encTopCRC.CRCPicAvgLeftBitsSum -= encTopCRC.CRCPicAvgLeftBitsSlot[encTopCRC.CRCPicAvgLeftBitsPrvIdx];
            encTopCRC.CRCPicAvgLeftBitsCnt--;
        }
        encTopCRC.CRCPicAvgLeftBitsSlot[encTopCRC.CRCPicAvgLeftBitsPrvIdx] = real_left_bits;
        encTopCRC.CRCPicAvgLeftBitsPrvIdx = (encTopCRC.CRCPicAvgLeftBitsPrvIdx + 1) % LEFTBITS_WIN;
    }

#if RESHAPE_BITS
#define GROUP_NUM 8
#define CHANGE_AMOUNT   0.4
    // closer to next I, more bits allocated
    const int change_num = (GopSize + GROUP_NUM - 1) / GROUP_NUM;
    const int change_idx = frameno / GROUP_NUM;
    double change_step = CHANGE_AMOUNT / (change_num - 1);
    /*if (change_idx == (change_num - 1))
    {
        const int last_frames_in_gop = GopSize - (change_num-1) * GROUP_NUM;
        change_step = change_step * GROUP_NUM / last_frames_in_gop;
    }*/
    const double rate_for_frame = 0.8 + change_idx * change_step;
    //ideal_bits = int(ideal_bits * rate_for_frame);
#else
    const double rate_for_frame = 1;
#endif

    // CBR
    // use smaller weight for GOPLeftBits, causu it's not good for steady target bits
    bits = int(pos_weight * target_bits + (1 - pos_weight) * gop_left_bits);
    printf("RC: P-Slice Bits(CP=%d) (%f*(ideal+VBuf/Gop), %f*gop_left) = (%d, %d, %d, %d)\n", encTopCRC.CRCPicCompMode,
           pos_weight, 1 - pos_weight, int(pos_weight * target_bits),
           int((1 - pos_weight) * gop_left_bits), VBuf_diff_bits / GopSize, old_gop_left);

    // reshape the GOP bits curve
    //bits = max(MIN_BITS_P, int(bits * rate_for_frame));
    bits = int(bits * rate_for_frame);
#endif

    return bits;
}

int find_lambda_idx(int Idx, const double estLambda_cur)
{
    // Find Optimal Idx according to estLambda
    double delta_more, delta_less;
    double estLambda_opt = encTopCRC.CRCLambdaTBL[Idx];
    if (estLambda_cur > estLambda_opt)
    {
        do
        {
            delta_more = estLambda_cur - estLambda_opt;
            Idx++;
            if (Idx > MAX_CRC_QP * LAMBDA_NUMQP)
            {
                Idx = MAX_CRC_QP * LAMBDA_NUMQP;
                break;
            }
            estLambda_opt = encTopCRC.CRCLambdaTBL[Idx];
        } while (estLambda_cur > estLambda_opt);

        delta_less = estLambda_opt - estLambda_cur;
        if (delta_more < delta_less)
            Idx--;
    }
    else    // (estLambda_cur < estLambda_opt)
    {
        do
        {
            delta_less = estLambda_opt - estLambda_cur;
            Idx--;
            if (Idx < MIN_CRC_QP * LAMBDA_NUMQP)
            {
                Idx = MIN_CRC_QP * LAMBDA_NUMQP;
                break;
            }
            estLambda_opt = encTopCRC.CRCLambdaTBL[Idx];
        } while (estLambda_cur < estLambda_opt);

        delta_more = estLambda_cur - estLambda_opt;
        if (delta_less < delta_more)
            Idx++;
    }
    return Idx;
}

int find_lambda_offset(int PicActualBits, int PicTargetBits, double enlarge_changes)
{
#define TBL_RANGE   (MAX_CRC_DQP_P*LAMBDA_NUMQP-1)
    //double lambda_tlb[TBL_RANGE] = { 1.1264, 1.2687, 1.4291, 1.6098, 1.8132, 2.0424, 2.3006, 2.5914, 2.9190, 3.2880, 3.7035, 4.1716 };
    double lambda_tlb[TBL_RANGE] = { 1.1264, 1.2687, 1.4291, 1.6098, 1.8132, 2.0424, 2.3006};
    double delta_more = 10, delta_less = 10000000;

    int sign, offset = 0, tbl_range;
    double ratio_c, conv_fac;
    if (PicActualBits > PicTargetBits)
    {
        ratio_c = (Double)PicActualBits / PicTargetBits;
        sign = 1;
        conv_fac = enlarge_changes * CONVERGE_FACTOR;
        tbl_range = 2 + encTopCRC.CRCLambdaIdxAccP;
        encTopCRC.CRCLambdaIdxAccP++;
        encTopCRC.CRCLambdaIdxAccN = VBR_EN ? 2 : 0;
    }
    else
    {
        ratio_c = (Double)PicTargetBits / PicActualBits;
        sign = -1;
        conv_fac = enlarge_changes * CONVERGE_FACTOR / 2;
        tbl_range = 1 + encTopCRC.CRCLambdaIdxAccN;
        encTopCRC.CRCLambdaIdxAccN++;
        encTopCRC.CRCLambdaIdxAccP = VBR_EN ? 2 : 0;
    }
    tbl_range = min(tbl_range, TBL_RANGE);
    ratio_c = (ratio_c - 1) * conv_fac + 1;

    // early terminate for small changes 
    if (ratio_c <= 1.015 && encTopCRC.CRCPicCompMode==0)
    {
        encTopCRC.CRCLambdaIdxAccN = VBR_EN ? 2 : 0;
        encTopCRC.CRCLambdaIdxAccP = VBR_EN ? 2 : 0;
        return 0;
    }

    if (ratio_c > lambda_tlb[tbl_range - 1])
    {
        offset = tbl_range - 1;
    }
    else
    {
        int idx;
        for (idx = 0; idx < tbl_range; idx++)
        {
            if (ratio_c <= lambda_tlb[idx])
            {
                delta_less = lambda_tlb[idx] - ratio_c;
                break;
            }
            else //if(ratio_c > lambda_tlb[idx])
            {
                delta_more = ratio_c - lambda_tlb[idx];
            }
        }
        offset = idx;
        if (delta_more < delta_less)
            offset--;
    }
    offset = offset + 1;
    printf("RC: Slice BPP ratio_c=%f, lambda_idx(old,new)=(%d,%d)\n", ratio_c, encTopCRC.CRCLambdaIdx, encTopCRC.CRCLambdaIdx+offset*sign);

#undef TBL_RANGE
    return offset*sign;
}

int find_sad_lambda_offset(int CurSADBits, int AvgSADBits)
{
#define TBL_RANGE   4
    //double lambda_tlb[TBL_RANGE] = { 1.1264, 1.2687, 1.4291, 1.6098, 1.8132, 2.0424, 2.3006, 2.5914, 2.9190, 3.2880, 3.7035, 4.1716 };
    double lambda_tlb[TBL_RANGE+1] = { 1.1264, 1.2687, 1.4291, 1.6098, 1.8132 };
    double delta_more = 10, delta_less = 10000000;

    int sign, offset = 0, tbl_range;
    double ratio_c, conv_fac;
    if (CurSADBits > AvgSADBits)
    {
        ratio_c = (Double)CurSADBits / AvgSADBits;
        sign = 1;
        conv_fac = 3*CONVERGE_FACTOR;
        tbl_range = TBL_RANGE - 1;

        ratio_c = (ratio_c - 1) * conv_fac + 1;
        // early terminate for small changes 
        if (ratio_c <= 1.075)
        {
            encTopCRC.CRCPicSADLambdaShift = 0;
            return 0;
        }
    }
    else
    {
        ratio_c = (Double)AvgSADBits / CurSADBits;
        sign = -1;
        conv_fac = 2*CONVERGE_FACTOR;
        tbl_range = TBL_RANGE - 1;

        ratio_c = (ratio_c - 1) * conv_fac + 1;
        // early terminate for small changes 
        if (ratio_c <= 1.05)
        {
            encTopCRC.CRCPicSADLambdaShift = -1;
            return -1;
        }
    }
    //tbl_range = min(tbl_range, TBL_RANGE);

    if (ratio_c > lambda_tlb[tbl_range - 1])
    {
        offset = tbl_range - 1;
    }
    else
    {
        int idx;
        for (idx = 0; idx < tbl_range; idx++)
        {
            if (ratio_c <= lambda_tlb[idx])
            {
                delta_less = lambda_tlb[idx] - ratio_c;
                break;
            }
            else //if(ratio_c > lambda_tlb[idx])
            {
                delta_more = ratio_c - lambda_tlb[idx];
            }
        }
        offset = idx;
        if (delta_more < delta_less)
            offset--;
    }
    offset = offset + 1;

    // limit range if trend not match with gop_left_bits
    if (encTopCRC.CRCPicCompMode == 1)  // lower bits
    {
        if (sign == -1)
            offset = min(1, offset);// max(0, offset - 1);
        else
            offset = Clip3(0, tbl_range, offset + 1);
    }
    else if (encTopCRC.CRCPicCompMode == 2) // raise bits
    {
        if (sign == 1)
            offset = min(1, offset);// max(0, offset - 1);
        else
            offset = Clip3(0, tbl_range, offset + 1);
    }

    printf("RC: Slice SAD ratio_c=%f, lambda_idx(old,new)=(%d,%d)\n", ratio_c, encTopCRC.CRCLambdaIdx, encTopCRC.CRCLambdaIdx + offset * sign);

#undef TBL_RANGE
    encTopCRC.CRCPicSADLambdaShift = offset * sign;
    return offset * sign;
}

int CRC_SAD_QPADJUST()
{
    const int frameno = encTopStat.m_encFrameNum;
    int lambdaOffset = 0;
    // Calculate SAD lambda offsets
    if (encTopCRC.CRCPicAvgSADCnt != 0 && frameno > SAD_WIN)
    {
        encTopCRC.CRCPicAvgSAD = encTopCRC.CRCPicAvgSADSum / encTopCRC.CRCPicAvgSADCnt;
        encTopCRC.CRCPicAvgLambdaIdx = encTopCRC.CRCPicAvgLambdaIdxSum / encTopCRC.CRCPicAvgSADCnt;
        encTopCRC.CRCPicAvgSADAlpha = encTopCRC.CRCPicAvgSADAlphaSum / encTopCRC.CRCPicAvgSADCnt;
        encTopCRC.CRCPicAvgSADBeta = encTopCRC.CRCPicAvgSADBetaSum / encTopCRC.CRCPicAvgSADCnt;
        encTopCRC.CRCPicAvgSADConst = encTopCRC.CRCPicAvgSADConstSum / encTopCRC.CRCPicAvgSADCnt;

        int AvgSADBitsA = int(encTopCRC.CRCPicAvgSADAlpha * encTopCRC.CRCPicAvgSAD / encTopCRC.CRCPicAvgLambdaIdx + encTopCRC.CRCPicAvgSADBeta);
        int AvgPrvCurSAD = int(0.75*encTopCRC.CRCPicSAD + 0.25*encTopCRC.CRCPicPrvSAD);
        double estLambdaIdxAvgA = (encTopCRC.CRCPicAvgSADAlpha * AvgPrvCurSAD) / (AvgSADBitsA - encTopCRC.CRCPicAvgSADBeta);
        int lambdaOffsetAvgA = (int(estLambdaIdxAvgA + 0.5) - encTopCRC.CRCLambdaIdx);
        int sign = (lambdaOffsetAvgA >= 0) ? 1 : -1;
        int val = abs(lambdaOffsetAvgA);

        if (val > 8)
        {
            //const int overshoot = val > 32;
            //val -= 8;
            if (sign > 0)    // large SAD
            {
                lambdaOffsetAvgA = sign * (val >> 3);
                //lambdaOffsetAvgA += 1;
                if ((encTopCRC.CRCPicCompMode == 2) || (encTopCRC.CRCPicPrvLambdaOffset < 0)) // raise bits
                    lambdaOffsetAvgA = max(0, lambdaOffsetAvgA - 2);   // SAD will raise bits for CRCPicCompMode=2 ?
                else if (encTopCRC.CRCPicPrvLambdaOffset > 0)
                    lambdaOffsetAvgA -= 1; // lower the effect
                    //lambdaOffsetAvgA = (lambdaOffsetAvgA + 1) >> 1; // lower the effect
            }
            else if (sign < 0)    // small SAD
            {
                lambdaOffsetAvgA = sign * (((val >> 3) + (val >> 4) + 1) >> 1);
                //lambdaOffsetAvgA -= 1;
                if ((encTopCRC.CRCPicCompMode == 1) || (encTopCRC.CRCPicPrvLambdaOffset > 0))// lower bits
                    lambdaOffsetAvgA = min(0, lambdaOffsetAvgA + 2);   // SAD will lower bits for CRCPicCompMode=1 ?
                else if (encTopCRC.CRCPicPrvLambdaOffset < 0)
                    lambdaOffsetAvgA += 1;   // lower the effect
                    //lambdaOffsetAvgA = lambdaOffsetAvgA >> 1;   // lower the effect
            }
        }
        else
        {
            lambdaOffsetAvgA = 0;
        }
        lambdaOffsetAvgA = Clip3(-3, 6, lambdaOffsetAvgA);
        printf("RC: SAD(Cur, PrvC, AvgA)=(%d, %d, %d) AvgASADBits=%d, lambda_idx(old, new, off)=(%d, %f, %d)\n", encTopCRC.CRCPicSAD, AvgPrvCurSAD,
               encTopCRC.CRCPicAvgSAD, AvgSADBitsA, encTopCRC.CRCLambdaIdx, estLambdaIdxAvgA, lambdaOffsetAvgA);
        lambdaOffset = lambdaOffsetAvgA;
    }

    // Avoid unreasonable values 
    int limit_SAD_min, limit_SAD_max, limit_SAD;
    limit_SAD = encTopCRC.CRCPicSAD;
    if (encTopCRC.CRCPicAvgSAD != 0)
    {
        limit_SAD_min = int(0.5 * encTopCRC.CRCPicAvgSAD);
        limit_SAD_max = int(1.5 * encTopCRC.CRCPicAvgSAD);
        limit_SAD = Clip3(limit_SAD_min, limit_SAD_max, limit_SAD);
    }

    // Sliding Window for Average SAD
    encTopCRC.CRCPicAvgSADSum += limit_SAD;
    encTopCRC.CRCPicAvgLambdaIdxSum += encTopCRC.CRCLambdaIdx;
    encTopCRC.CRCPicAvgSADAlphaSum += encTopCRC.CRCPicSADAlpha;
    encTopCRC.CRCPicAvgSADBetaSum += encTopCRC.CRCPicSADBeta;
    encTopCRC.CRCPicAvgSADConstSum += encTopCRC.CRCPicSADConst;
    encTopCRC.CRCPicAvgSADCnt++;
    if (encTopCRC.CRCPicAvgSADCnt > SAD_WIN)
    {
        encTopCRC.CRCPicAvgSADSum -= encTopCRC.CRCPicAvgSADSlot[encTopCRC.CRCPicAvgSADPrvIdx];
        encTopCRC.CRCPicAvgLambdaIdxSum -= encTopCRC.CRCPicAvgLambdaIdxSlot[encTopCRC.CRCPicAvgSADPrvIdx];
        encTopCRC.CRCPicAvgSADAlphaSum -= encTopCRC.CRCPicAvgSADAlphaSlot[encTopCRC.CRCPicAvgSADPrvIdx];
        encTopCRC.CRCPicAvgSADBetaSum -= encTopCRC.CRCPicAvgSADBetaSlot[encTopCRC.CRCPicAvgSADPrvIdx];
        encTopCRC.CRCPicAvgSADConstSum -= encTopCRC.CRCPicAvgSADConstSlot[encTopCRC.CRCPicAvgSADPrvIdx];
        encTopCRC.CRCPicAvgSADCnt--;
    }
    encTopCRC.CRCPicAvgSADSlot[encTopCRC.CRCPicAvgSADPrvIdx] = limit_SAD;
    encTopCRC.CRCPicAvgLambdaIdxSlot[encTopCRC.CRCPicAvgSADPrvIdx] = encTopCRC.CRCLambdaIdx;
    encTopCRC.CRCPicAvgSADAlphaSlot[encTopCRC.CRCPicAvgSADPrvIdx] = encTopCRC.CRCPicSADAlpha;
    encTopCRC.CRCPicAvgSADBetaSlot[encTopCRC.CRCPicAvgSADPrvIdx] = encTopCRC.CRCPicSADBeta;
    encTopCRC.CRCPicAvgSADConstSlot[encTopCRC.CRCPicAvgSADPrvIdx] = encTopCRC.CRCPicSADConst;
    encTopCRC.CRCPicAvgSADPrvIdx = (encTopCRC.CRCPicAvgSADPrvIdx + 1) % SAD_WIN;

    return lambdaOffset;
}

int Check_CPB_Fullness(int estBits)
{
    const int frameno = encTopStat.m_encFrameNum;
    const int GopSize = encTopCRC.CRCGopSize;
    int newBits = estBits;
    //if (encTopCRC.CRCCpbEn)
    {
        // Find Lower Bound
#if VBR_EN
        int LowerBound = encTopCRC.CRCBufferingRate;
#else
        int LowerBound = encTopCRC.CRCBufferingRate;
        /*const int next_frm_pos = (frameno + 1) % GopSize;
        const int left_frm = GopSize - next_frm_pos;
        int GopBits = encTopCRC.CRCGopLeftBits - estBits;
        if (next_frm_pos == 0)
        {
            GopBits += CRC_Init_GOP();  // this is not very accurate ?
        }
        int LowerBound = (GopBits + left_frm - 1) / left_frm;
        LowerBound = int(LowerBound * 0.2 + encTopCRC.CRCTargetAvgBit * 0.8);*/
#endif
        const int CPBRecover = (int)(encTopCRC.CRCCpbSize * 0.5 + 0.5);

        // Control CPB fullness
        /*if (encTopCRC.CRCCpbEn)
        {
            //if (isISlice && (encTopCRC.CRCCpbState < int(encTopCRC.CRCCpbSize * 0.25)))
            if (isISlice && (encTopCRC.CRCCpbState < int(encTopCRC.CRCCpbSize * 0.75)))
            {
                encTopCRC.CRCPicEstHADIQP += 1; // lower actual bits
                const int HAD = encTopStat.FrameStat[frameno].m_avgHAD;
                const double alpha = encTopStat.HAD_ALPHA[encTopCRC.CRCPicEstHADIQP];
                const double beta = encTopStat.HAD_BETA[encTopCRC.CRCPicEstHADIQP];
                const int bits = int(alpha * HAD + beta);
                const int max_Ibits = encTopCRC.CRCGopTargetBits - int((GopSize - 1) * MIN_BITS_P * encTopCRC.CRCHadScale);
                const int min_Ibits = int(MIN_BITS_I * encTopCRC.CRCHadScale);
                encTopCRC.CRCPicEstHADIBits = estBits = Clip3(min_Ibits, max_Ibits, bits);
                printf("RC: CPB State = %3.1f%%, QP=%d\n", (float)encTopCRC.CRCCpbState * 100 / encTopCRC.CRCCpbSize, encTopCRC.CRCPicEstHADIQP);
            }
        }*/

        int estimatedCpbFullness = encTopCRC.CRCCpbState + encTopCRC.CRCBufferingRate;
        // prevent overflow (estBits too low)
        if ((estimatedCpbFullness - estBits) > int(encTopCRC.CRCCpbSize * 0.9))
        {
            if (encTopCRC.CRCCpbEn)
            {
                newBits = estimatedCpbFullness - int(encTopCRC.CRCCpbSize * 0.9);
                //newBits = min(int(1.5 * old_bits), estimatedCpbFullness - int(encTopCRC.CRCCpbSize * 0.9));
                printf("RC: CPB overflow estBits from %d to %d\n", estBits, newBits);
            }
            else
            {
                encTopStat.FrameStat[frameno].cpb_overflow = 1;
                encTopCRC.CRCCpbState = CPBRecover;
            }
        }

        estimatedCpbFullness -= encTopCRC.CRCBufferingRate;
        // prevent underflow (estBits too high)
        if ((estimatedCpbFullness - estBits) < LowerBound)
        {
            if (encTopCRC.CRCCpbEn)
            {
                newBits = max(int(MIN_BITS_P * encTopCRC.CRCHadScale), estimatedCpbFullness - LowerBound);
                //newBits = max(int(0.5 * old_bits), estimatedCpbFullness - LowerBound);
                printf("RC: CPB underflow estBits from %d to %d\n", estBits, newBits);
            }
            else
            {
                encTopStat.FrameStat[frameno].cpb_overflow = -1;
                encTopCRC.CRCCpbState = LowerBound;
            }
        }
    }
    return newBits;
}

void CRC_Update_LambdaIdx(int PicTargetBits, int PicActualBits)
{
    const int frameno = encTopStat.m_encFrameNum;
    const int GopSize = encTopCRC.CRCGopSize;
    const int frm_pos = frameno % GopSize;
    const int left_frm = GopSize - frm_pos;
    int is_update = 1;

    // Control CPB fullness
    if (encTopCRC.CRCCpbEn)
    {
        int frame_idx = left_frm >> 3;
        const double rate_pos_h[7] = { 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.4 };
        const double rate_pos_l[7] = { 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.2 };
        if ((encTopCRC.CRCCpbState > int(encTopCRC.CRCCpbSize * rate_pos_h[frame_idx])))  // CPB increase too fast = actual bits too low
        {
            // applied only when gop_left_bits trend matches
            if (encTopCRC.CRCPicCompMode != 1)  // 0 or 2
            {
                encTopCRC.CRCPicCompMode = 2;   // raise actual bits
                PicTargetBits = max(PicTargetBits, encTopCRC.CRCBufferingRate);
                encTopCRC.CRCLambdaIdxAccN += 1; // (encTopCRC.CRCPicCompMode == 2);
                //encTopCRC.CRCLambdaIdx -= 1;
                printf("RC: CPB State = %3.1f%%, Target = %3.1f%%, CP_R=%d\n", (float)encTopCRC.CRCCpbState * 100 / encTopCRC.CRCCpbSize, 100 * rate_pos_h[frame_idx], encTopCRC.CRCPicCompMode);
            }
        }
        else if ((encTopCRC.CRCCpbState < int(encTopCRC.CRCCpbSize * rate_pos_l[frame_idx])))  // CPB increase too slow = actual bits too high
        {
            // applied only when gop_left_bits trend matches
            if (encTopCRC.CRCPicCompMode != 2)  // 0 or 1
            {
                encTopCRC.CRCPicCompMode = 1;   // lower actual bits
                PicTargetBits = min(PicTargetBits, encTopCRC.CRCBufferingRate);
                encTopCRC.CRCLambdaIdxAccP += 1;
                //encTopCRC.CRCLambdaIdx += 1;
                printf("RC: CPB State = %3.1f%%, Target = %3.1f%%, CP_L=%d\n", (float)encTopCRC.CRCCpbState * 100 / encTopCRC.CRCCpbSize, 100 * rate_pos_l[frame_idx], encTopCRC.CRCPicCompMode);
            }
        }
    }

    // keep compress more/less, do not update parameters
    Double enlarge_changes = VBR_EN ? 1 : 1;
    if (encTopCRC.CRCPicCompMode != 0)
    {
        if (encTopCRC.CRCPicCompMode == 1 && (PicTargetBits > PicActualBits))   // gop_left_bits too low, lower actual bits
        {
            encTopCRC.CRCLambdaIdxAccP = VBR_EN ? 2 : 0;
            encTopCRC.CRCLambdaIdxAccN = VBR_EN ? 2 : 0;
            is_update = 0;
#if CBR_SAD_EN
            encTopCRC.CRCPicPrvLambdaOffset = 1;    // for record CRCPicCompMode=1
#endif
        }
        else if (encTopCRC.CRCPicCompMode == 2 && (PicTargetBits < PicActualBits))  // gop_left_bits too high, raise actual bits
        {
            encTopCRC.CRCLambdaIdxAccP = VBR_EN ? 2 : 0;
            encTopCRC.CRCLambdaIdxAccN = VBR_EN ? 2 : 0;
            is_update = 0;
#if CBR_SAD_EN
            encTopCRC.CRCPicPrvLambdaOffset = -1;   // for record CRCPicCompMode=2
#endif
        }
        else    // enlarge the changes ?
            enlarge_changes *= 2;
    }

    // Ctu Average offset
    //encTopCRC.CRCLambdaIdx += encTopCRC.CRCCtuLambdaIdxOffset;

    // Find Lambda Idx offset
    if (is_update)
    {
        const int lambda_offset = find_lambda_offset(PicActualBits, PicTargetBits, enlarge_changes);
        encTopCRC.CRCLambdaIdx += lambda_offset;
#if CBR_SAD_EN
        encTopCRC.CRCPicPrvLambdaOffset = lambda_offset;
#endif
    }
}

void CRC_Init_PIC(int isISlice)
{
    const int frameno = encTopStat.m_encFrameNum;
    const int GopSize = encTopCRC.CRCGopSize;
    int estBits;

    // Derek debug
    if (frameno == 44)
        estBits = 0;

    if (isISlice)    // I-Slice
    {
#if !VBR_EN
        encTopCRC.CRCGopTargetBits = encTopCRC.CRCGopLeftBits = CRC_Init_GOP();

        // reset Sliding Window for Average LeftBits
        encTopCRC.CRCPicAvgLeftBitsCnt = 0;
        encTopCRC.CRCPicAvgLeftBitsPrvIdx = 0;
        encTopCRC.CRCPicAvgLeftBitsSum = 0;
#endif
        estBits = CRC_Estimate_IBits(encTopCRC.CRCPrevPQP);
        encTopCRC.CRCPicTargetBits = Check_CPB_Fullness(estBits);
        // re-calculate IQP if IBits is cut
        if (encTopCRC.CRCPicTargetBits != estBits)
        {
            encTopCRC.CRCPicEstHADIQP = CRC_Estimate_IBits2(encTopCRC.CRCPicEstHADIQP, encTopCRC.CRCPicTargetBits);
            //encTopCRC.CRCPicEstHADIQP = Clip3(encTopCRC.CRCPrevPQP - MAX_CRC_DQP_I_N - 1, encTopCRC.CRCPrevPQP + MAX_CRC_DQP_I_P + 1, encTopCRC.CRCPicEstHADIQP);   // prevent IP flickering
            //encTopCRC.CRCPicEstHADIQP = Clip3(MIN_CRC_QP, MAX_CRC_QP, encTopCRC.CRCPicEstHADIQP);
            const int HAD = encTopStat.FrameStat[frameno].m_avgHAD;
            const double alpha = encTopStat.HAD_ALPHA[encTopCRC.CRCPicEstHADIQP];
            const double beta = encTopStat.HAD_BETA[encTopCRC.CRCPicEstHADIQP];
            const int bits = int(alpha * HAD + beta);
#if VBR_EN
            const int max_Ibits = int(MAX_BITS_I * encTopCRC.CRCHadScale);
#else
            const int max_Ibits = encTopCRC.CRCGopTargetBits - int((GopSize - 1) * MIN_BITS_P * encTopCRC.CRCHadScale);
#endif
            const int min_Ibits = int(MIN_BITS_I * encTopCRC.CRCHadScale);
            encTopCRC.CRCPicTargetBits = Clip3(min_Ibits, max_Ibits, bits);
        }
        encTopCRC.CRCPicEstHADIBits = encTopCRC.CRCPicTargetBits;

        // Find Lamdba
        const Int sliceQP = encTopCRC.CRCPicQP = encTopCRC.CRCPicEstHADIQP;
        encTopCRC.CRCLambdaIdx = sliceQP * LAMBDA_NUMQP;
        encTopCRC.CRCPicLambda = encTopCRC.CRCLambdaTBL[encTopCRC.CRCLambdaIdx];
    }
    else
    {
#if VBR_EN
        encTopCRC.CRCPicEstActBits = encTopStat.FrameStat[frameno - 1].frameBits;
#endif
        estBits = CRC_Estimate_PBits();
        encTopCRC.CRCPicTargetBits = Check_CPB_Fullness(estBits);
         
        // Find Lamdba
        if (frameno == 1)   // fisrt P frame
        {
            // Derek, changes with QP?
            Double estLambda = I_TO_P_EST * encTopCRC.CRCPicLambda;
            encTopCRC.CRCLambdaIdx = find_lambda_idx(encTopCRC.CRCLambdaIdx, estLambda);
        }
        else if ((frameno % encTopCRC.CRCGopSize) == 1)
        {
            const Double bpp = encTopCRC.CRCPicTargetBits / (Double)encTopCRC.CRCPicPixels;
            Double estLambda = sqrt(encTopCRC.CRCPicConstPAvg / bpp);
            encTopCRC.CRCLambdaIdx = find_lambda_idx(encTopCRC.CRCLambdaIdx, estLambda);
        }
        else
        {
#if VBR_EN
            const int PrvActualBits = encTopCRC.CRCPicEstActBits;
#else
            const int PrvActualBits = encTopStat.FrameStat[frameno - 1].frameBits;
#endif
            CRC_Update_LambdaIdx(encTopCRC.CRCPicTargetBits, PrvActualBits);
        }

#if CBR_SAD_EN
        encTopCRC.CRCPicSAD = encTopStat.FrameStat[frameno].m_avgSAD;
        if (encTopCRC.CRCPicSAD == 0)
        {
            int lambda_offset = -1;
            encTopCRC.CRCPicSADLambdaShift = lambda_offset;
            encTopCRC.CRCLambdaIdx += lambda_offset;
        }
        else
        {
            int lambda_offset = CRC_SAD_QPADJUST();
            encTopCRC.CRCPicSADLambdaShift = lambda_offset;
            encTopCRC.CRCLambdaIdx += lambda_offset;
        }
        encTopCRC.CRCPicPrvLambdaOffset = 0;
#endif

        // Find QP
        int estQP = encTopCRC.CRCLambdaIdx / LAMBDA_NUMQP;
        encTopCRC.CRCPicLambda = encTopCRC.CRCLambdaTBL[encTopCRC.CRCLambdaIdx];
        encTopCRC.CRCPicQP = Clip3(encTopCRC.CRCPrevPQP - MAX_CRC_DQP_P, encTopCRC.CRCPrevPQP + MAX_CRC_DQP_P, estQP);
        encTopCRC.CRCPicQP = Clip3(MIN_CRC_QP, MAX_CRC_QP, encTopCRC.CRCPicQP);
#if VBR_EN
        encTopCRC.CRCPicQP = Clip3(encTopCRC.CRCPicVBRMinQP, encTopCRC.CRCPicVBRMaxQP, encTopCRC.CRCPicQP);
#endif
        // recalculate if QP is clipped
        if (estQP != encTopCRC.CRCPicQP)
        {
            double old_lambda = encTopCRC.CRCPicLambda;
            encTopCRC.CRCLambdaIdx = (estQP > encTopCRC.CRCPicQP) ? ((encTopCRC.CRCPicQP + 1) * LAMBDA_NUMQP - 1) : (encTopCRC.CRCPicQP * LAMBDA_NUMQP);
            encTopCRC.CRCPicLambda = encTopCRC.CRCLambdaTBL[encTopCRC.CRCLambdaIdx];
            int new_bits = int(encTopCRC.CRCPicTargetBits * old_lambda * old_lambda / encTopCRC.CRCPicLambda / encTopCRC.CRCPicLambda);
            printf("RC: QP clipped from %d to %d, estBits from %d to %d\n", estQP, encTopCRC.CRCPicQP, encTopCRC.CRCPicTargetBits, new_bits);
        }
    }
    encTopCRC.CRCCtuLambdaIdx = encTopCRC.CRCLambdaIdx;
    encTopCRC.CRCCtuLambda = encTopCRC.CRCPicLambda;
    encTopCRC.CRCCtuQP = encTopCRC.CRCPicQP;
    encTopCRC.CRCCtuAccumBits = 0;
}

void CRC_Update_Params(int isISlice)
{
    const int frameno = encTopStat.m_encFrameNum;
    const int PicActualBits = encTopStat.FrameStat[frameno].frameBits;

    // Update Ctu Average LambdaIdx QP Lambda
    if (encTopCRC.CRCCtuEn)
    {
        int oldLIdx = encTopCRC.CRCLambdaIdx;
        double oldLambda = encTopCRC.CRCPicLambda;
        int oldQP = encTopCRC.CRCPicQP;
        
        int newIdx = encTopCRC.CRCCtuAvgLambdaIdx / encTopCRC.CRCCtuValidCnt;
        double newLambda = encTopCRC.CRCLambdaTBL[newIdx];
        int newQP = newIdx / LAMBDA_NUMQP;

        //if (newIdx != oldLIdx)
        //    encTopCRC.CRCCtuLambdaScale = pow(1.0613235, newIdx - oldLIdx);

        printf("RC: CTU Post Update QP(old,new)=(%d, %d), LIdx=(%d, %d), Lambda=(%f, %f), Scale=%f\n", oldQP, newQP, oldLIdx, newIdx, oldLambda, newLambda, encTopCRC.CRCCtuLambdaScale);

        //encTopCRC.CRCLambdaIdx = newIdx;
        //encTopCRC.CRCPicLambda = newLambda;
        //encTopCRC.CRCPicQP = newQP;
        encTopCRC.CRCCtuLambdaIdxOffset = encTopCRC.CRCLambdaIdx - newIdx;

        encTopCRC.CRCCtuAvgLambdaIdx = 0;
        //encTopCRC.CRCCtuAvgQP = 0;
        encTopCRC.CRCCtuValidCnt = 0;
    }

    // Update I-slice HAD params
    if (isISlice)
    {
        const int QP = encTopCRC.CRCPicQP;  // ==encTopCRC.CRCPicEstHADIQP
        const double alpha = encTopStat.HAD_ALPHA[QP];
        const double beta = encTopStat.HAD_BETA[QP];
        const int delta_bits = PicActualBits - encTopCRC.CRCPicEstHADIBits;
        const int HAD = encTopStat.FrameStat[frameno].m_avgHAD;
        const double new_beta = beta + 0.5 * delta_bits;
        const double new_alpha = alpha + (0.5 * delta_bits) / HAD;
        encTopStat.HAD_ALPHA[QP] = new_alpha;
        encTopStat.HAD_BETA[QP] = new_beta;
        printf("RC: I-Slice HAD(QP%d) alpha(old,new)=(%f,%f), beta(old,new)=(%f,%f)\n", QP, alpha, new_alpha, beta, new_beta);

        encTopStat.FrameStat[frameno].pic_HAD_alpha = new_alpha;
        encTopStat.FrameStat[frameno].pic_HAD_beta = new_beta;

        // Derek, Update other QPs to better guess next GOP-I ?
        if (frameno == 0)
        {
            // HAD alpha is proximately linear with QP ?, beta is not linear with QP
            const double alpha_rate = new_alpha / alpha;
            const double beta_rate = (beta==0) ? alpha_rate : (new_beta / beta);
            for (int i = MIN_CRC_QP; i < QP; i++)
            {
                encTopStat.HAD_ALPHA[i] *= alpha_rate;
                encTopStat.HAD_BETA[i] *= beta_rate;
            }
            for (int i = QP+1; i <= MAX_CRC_QP; i++)
            {
                encTopStat.HAD_ALPHA[i] *= alpha_rate;
                encTopStat.HAD_BETA[i] *= beta_rate;
            }

            // reset Sliding Window for Previous IBits
            encTopCRC.CRCPicAvgICnt = 0;
            encTopCRC.CRCPicAvgIPrvIdx = 0;
            encTopCRC.CRCPicAvgIBitsSum = 0;
        }

        // Sliding Window for Previous IBits
        encTopCRC.CRCPicAvgIBitsSum += PicActualBits;
        encTopCRC.CRCPicAvgICnt++;
        if (encTopCRC.CRCPicAvgICnt > SMOOTH_I_WIN)
        {
            encTopCRC.CRCPicAvgIBitsSum -= encTopCRC.CRCPicAvgIBitsSlot[encTopCRC.CRCPicAvgIPrvIdx];
            encTopCRC.CRCPicAvgICnt--;
        }
        encTopCRC.CRCPicAvgIBitsSlot[encTopCRC.CRCPicAvgIPrvIdx] = PicActualBits;
        encTopCRC.CRCPicAvgIPrvIdx = (encTopCRC.CRCPicAvgIPrvIdx + 1) % SMOOTH_I_WIN;

        encTopCRC.CRCPicAvgIBits = encTopCRC.CRCPicAvgIBitsSum / encTopCRC.CRCPicAvgICnt;

#if VBR_EN
        encTopCRC.CRCPicPrvIBits = PicActualBits;
#else
        encTopCRC.CRCPicIdealBits2 = encTopCRC.CRCPicIdealBits - (PicActualBits - encTopCRC.CRCTargetAvgBit) / (encTopCRC.CRCGopSize - 1);  // g_RCSmoothWindowSize
#endif
    }
#if SAD_EN
    else if(encTopCRC.CRCPicSAD != 0)//if (!isISlice)
    {
        const int lambdaIdx = encTopCRC.CRCLambdaIdx;
        const int SAD = encTopCRC.CRCPicSAD;
        const double alpha = encTopCRC.CRCPicSADAlpha;
        const double beta = encTopCRC.CRCPicSADBeta;
        //const int cnst = encTopCRC.CRCPicSADConst;
        const int est_bits = int(alpha * SAD / lambdaIdx + beta);
        //const int est_bits = int(alpha * SAD / lambdaIdx + beta * SAD + cnst);
        //const int est_bits = int(beta * SAD - alpha * SAD * lambdaIdx + cnst);
        const int delta_bits = PicActualBits - est_bits;

        encTopCRC.CRCPicSADAlpha = alpha + (double)delta_bits * lambdaIdx / SAD / 2;
        encTopCRC.CRCPicSADBeta = beta + (double)delta_bits / 2;
        //encTopCRC.CRCPicSADAlpha = alpha + (double)delta_bits * lambdaIdx / SAD / 3;
        //encTopCRC.CRCPicSADBeta = beta + (double)delta_bits / SAD / 3;
        //encTopCRC.CRCPicSADConst = cnst + (delta_bits / 3);

        // keep parameters stable
        if (encTopCRC.CRCPicSADBeta < 0)
        {
            encTopCRC.CRCPicSADAlpha += encTopCRC.CRCPicSADBeta * lambdaIdx / SAD;
            encTopCRC.CRCPicSADBeta = 0;
        }
        if (encTopCRC.CRCPicSADAlpha < 0)
        {
            const double alpha_add = 500;
            encTopCRC.CRCPicSADBeta += (encTopCRC.CRCPicSADAlpha - alpha_add) * SAD / lambdaIdx;
            encTopCRC.CRCPicSADAlpha = alpha_add;
        }
        
        /*// keep alpha > 0
        double alpha_delta = (double)delta_bits / lambdaIdx / SAD / 3;
        if (alpha_delta > alpha)
        {
            double diff = alpha_delta - alpha + 0.5;
            encTopCRC.CRCPicSADAlpha = 0.5;
            encTopCRC.CRCPicSADBeta = beta + ((double)delta_bits / SAD / 3) + (diff * lambdaIdx / 2);
            encTopCRC.CRCPicSADConst = cnst + (delta_bits / 3) + (diff * lambdaIdx * SAD / 2);
        }
        else
        {
            encTopCRC.CRCPicSADAlpha = alpha - (double)delta_bits / lambdaIdx / SAD / 3;
            encTopCRC.CRCPicSADBeta = beta + (double)delta_bits / SAD / 3;
            encTopCRC.CRCPicSADConst = cnst + (delta_bits / 3);
        }*/

        encTopCRC.CRCPicPrvSAD = SAD;
        //encTopCRC.CRCPicPrvLambdaIdx = encTopCRC.CRCLambdaIdx;

        int new_bits = int(encTopCRC.CRCPicSADAlpha * SAD / lambdaIdx + encTopCRC.CRCPicSADBeta);
        printf("RC: P-Slice SAD_%d(LIdx%d) alpha(old,new)=(%f,%f), beta(old,new)=(%f,%f), const(old,new)=(%d,%d), new_bits=%d\n", SAD, lambdaIdx, alpha, encTopCRC.CRCPicSADAlpha,
               beta, encTopCRC.CRCPicSADBeta, /*cnst*/0, encTopCRC.CRCPicSADConst, new_bits);
    }
#endif

    // Update GOP Params
    encTopCRC.CRCPrevPQP = encTopCRC.CRCPicQP;
    encTopCRC.CRCFramesLeft--;
    encTopCRC.CRCGopDeltaBits += (PicActualBits - encTopCRC.CRCTargetAvgBit);
    if (isISlice)
    {
        encTopCRC.CRCGopTargetVBuf = encTopCRC.CRCGopDeltaBits;
        encTopCRC.CRCGopAvgVBuf = encTopCRC.CRCGopDeltaBits / encTopCRC.CRCGopSize;
        //encTopCRC.CRCPicIdealBits2 = encTopCRC.CRCTargetAvgBit - encTopCRC.CRCGopAvgVBuf;
    }
    else
    {
        encTopCRC.CRCGopTargetVBuf -= encTopCRC.CRCGopAvgVBuf;
    }
    encTopCRC.CRCGopLeftBits -= PicActualBits;  // Derek, for intra frame, HM uses estimated bits to update ?

    // update CPB state
    encTopCRC.CRCCpbState += encTopCRC.CRCBufferingRate - PicActualBits;

    // record status
    encTopStat.FrameStat[frameno].cpb_state = encTopCRC.CRCCpbState;
    encTopStat.FrameStat[frameno].cpb_vbuf = encTopCRC.CRCGopDeltaBits;
    encTopStat.FrameStat[frameno].est_vbuf = encTopCRC.CRCGopTargetVBuf;
    encTopStat.FrameStat[frameno].pic_lambda = encTopCRC.CRCPicLambda;
    encTopStat.FrameStat[frameno].target_bits = encTopCRC.CRCPicTargetBits;
    encTopStat.FrameStat[frameno].gop_left_bits = encTopCRC.CRCGopLeftBits;

    if (isISlice)
    {
        encTopCRC.CRCLambdaIdxAccP = VBR_EN ? 2 : 0;
        encTopCRC.CRCLambdaIdxAccN = VBR_EN ? 2 : 0;
    }
    else
    {
#if CBR_SAD_EN
        // restore lambda_offset
        int sign = (encTopCRC.CRCPicSADLambdaShift > 0) ? 1 : -1;
        int val = abs(encTopCRC.CRCPicSADLambdaShift);
        val = sign * (val >> 1);
        encTopCRC.CRCLambdaIdx -= val;
        //encTopCRC.CRCLambdaIdx -= encTopCRC.CRCPicSADLambdaShift;
#endif

        // Sliding Window for ConstP
        encTopCRC.CRCPicConstP = encTopCRC.CRCPicLambda * encTopCRC.CRCPicLambda * PicActualBits / encTopCRC.CRCPicPixels;
        printf("RC: P-Slice ConstP = %f\n", encTopCRC.CRCPicConstP);    // watch if it breaks
        encTopCRC.CRCPicConstPSum += encTopCRC.CRCPicConstP;
        encTopCRC.CRCPicConstPBitsSum += PicActualBits;
        encTopCRC.CRCPicConstPCnt++;
        if (encTopCRC.CRCPicConstPCnt > CONSTP_WIN)
        {
            encTopCRC.CRCPicConstPSum -= encTopCRC.CRCPicConstPSlot[encTopCRC.CRCPicConstPPrvIdx];
            encTopCRC.CRCPicConstPBitsSum -= encTopCRC.CRCPicConstPBitsSlot[encTopCRC.CRCPicConstPPrvIdx];
            encTopCRC.CRCPicConstPCnt--;
        }
        encTopCRC.CRCPicConstPSlot[encTopCRC.CRCPicConstPPrvIdx] = encTopCRC.CRCPicConstP;
        encTopCRC.CRCPicConstPBitsSlot[encTopCRC.CRCPicConstPPrvIdx] = PicActualBits;
        encTopCRC.CRCPicConstPPrvIdx = (encTopCRC.CRCPicConstPPrvIdx + 1) % CONSTP_WIN;

#if VBR_EN
        encTopCRC.CRCPicConstPBitsAvg = encTopCRC.CRCPicConstPBitsSum / encTopCRC.CRCPicConstPCnt;
#endif

        // Derek, Update next I parameters ?
        const int is_lastGOP = ((frameno + 1) % encTopCRC.CRCGopSize) == 0;
        if (is_lastGOP)
        {
            encTopCRC.CRCPicConstPAvg = encTopCRC.CRCPicConstPSum / encTopCRC.CRCPicConstPCnt;
            encTopCRC.CRCPicConstPBitsAvg = encTopCRC.CRCPicConstPBitsSum / encTopCRC.CRCPicConstPCnt;
            encTopCRC.CRCPicConstI = P_TO_I_EST * encTopCRC.CRCPicConstPAvg * encTopCRC.CRCPicAvgIBits / encTopCRC.CRCPicConstPBitsAvg; // encTopCRC.CRCPicIdealBits2
            printf("RC: P-Slice (ConstI, AvgConstP) = (%f, %f)\n", encTopCRC.CRCPicConstI, encTopCRC.CRCPicConstPAvg);
        }
    }
}


#if 0
int find_ctu_lambda_offset(int PicActualBits, int PicTargetBits)
{
    //double lambda_tlb[TBL_RANGE] = { 1.1264, 1.2687, 1.4291, 1.6098, 1.8132, 2.0424, 2.3006, 2.5914, 2.9190, 3.2880, 3.7035, 4.1716 };
    double lambda_tlb[7] = { 1.1264, 1.2687, 1.4291, 1.6098, 1.8132, 2.0424, 2.3006 };
    double delta_more = 10, delta_less = 10000000;

    int sign, offset = 0, tbl_range;
    double ratio_c, conv_fac;
    if (PicActualBits > PicTargetBits)
    {
        ratio_c = (Double)PicActualBits / PicTargetBits;
        sign = 1;
        conv_fac = CONVERGE_FACTOR;
        tbl_range = 2;
    }
    else
    {
        ratio_c = (Double)PicTargetBits / PicActualBits;
        sign = -1;
        conv_fac = CONVERGE_FACTOR / 2;
        tbl_range = 1;
    }
    ratio_c = (ratio_c - 1) * conv_fac + 1;

    // early terminate for small changes 
    if (ratio_c <= 1.015 && encTopCRC.CRCCtuCompMode == 0)
    {
        return 0;
    }

    if (ratio_c > lambda_tlb[tbl_range - 1])
    {
        offset = tbl_range - 1;
    }
    else
    {
        int idx;
        for (idx = 0; idx < tbl_range; idx++)
        {
            if (ratio_c <= lambda_tlb[idx])
            {
                delta_less = lambda_tlb[idx] - ratio_c;
                break;
            }
            else //if(ratio_c > lambda_tlb[idx])
            {
                delta_more = ratio_c - lambda_tlb[idx];
            }
        }
        offset = idx;
        if (delta_more < delta_less)
            offset--;
    }
    offset = offset + 1;

    return offset * sign;
}

void CRC_Estimate_CTU_Lambda(TComPic* pcPic, int totalCtuInSlice, int CtuLeftNum)
{
    int totalCtuNum = pcPic->getNumberOfCtusInFrame();
    int TargetSliceBits = encTopCRC.CRCPicTargetBits * totalCtuInSlice / totalCtuNum; // Derek, replace TargetBits with real left bits ?
    int SliceLeftBits = TargetSliceBits - encTopCRC.CRCCtuAccumBits;
    int avgCtuBits = encTopCRC.CRCPicTargetBits / totalCtuNum;
    int avgLeftBits = SliceLeftBits / CtuLeftNum;

    // detect if compress more
    int old_ctu_left = avgLeftBits;
    if (avgLeftBits < int((1 - COMPRESS_FACTOR) * avgCtuBits))
    {
        encTopCRC.CRCCtuCompMode = 1;
        // clip avgLeftBits if unreasonable
        if (avgLeftBits < int(0.8 * avgCtuBits))
            avgLeftBits = int(0.8 * avgCtuBits);
    }
    else if (avgLeftBits > int((1 + COMPRESS_FACTOR) * avgCtuBits))
    {
        encTopCRC.CRCCtuCompMode = 2;
        // clip gop_left_bits if unreasonable
        if (avgLeftBits > int(1.2 * avgCtuBits))
            avgLeftBits = int(1.2 * avgCtuBits);
    }
    else
        encTopCRC.CRCCtuCompMode = 0;

    encTopCRC.CRCCtuTargetBits = (avgCtuBits + avgLeftBits) >> 1;

    // keep compress more/less, do not update parameters
    int is_update = 1;
    int PicTargetBits = encTopCRC.CRCCtuTargetBits;
    int PicActualBits = (encTopCRC.CRCCtuActualBits==0) ? PicTargetBits : encTopCRC.CRCCtuActualBits; // previous ctu
    if (encTopCRC.CRCCtuCompMode != 0)
    {
        if (encTopCRC.CRCCtuCompMode == 1 && (PicTargetBits > PicActualBits))   // left bits too low, lower actual bits
        {
            is_update = 0;
        }
        else if (encTopCRC.CRCCtuCompMode == 2 && (PicTargetBits < PicActualBits))  // left bits too high, raise actual bits
        {
            is_update = 0;
        }
    }

    // Find Lambda Idx offset
    if (is_update)
    {
        const int lambda_offset = find_ctu_lambda_offset(PicActualBits, PicTargetBits);
        encTopCRC.CRCCtuLambdaIdx += lambda_offset;

        int max_idx = encTopCRC.CRCLambdaIdx + 8;
        int min_idx = encTopCRC.CRCLambdaIdx - 8;
        encTopCRC.CRCCtuLambdaIdx = Clip3(min_idx, max_idx, encTopCRC.CRCCtuLambdaIdx);
        encTopCRC.CRCCtuLambda = encTopCRC.CRCLambdaTBL[encTopCRC.CRCCtuLambdaIdx];

        int estQP = encTopCRC.CRCCtuLambdaIdx / LAMBDA_NUMQP;
        encTopCRC.CRCCtuQP = Clip3(encTopCRC.CRCPicQP - MAX_CRC_DQP_P, encTopCRC.CRCPicQP + MAX_CRC_DQP_P, estQP);
        encTopCRC.CRCCtuQP = Clip3(MIN_CRC_QP, MAX_CRC_QP, encTopCRC.CRCCtuQP);
    }
}
#endif