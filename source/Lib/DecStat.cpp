

#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "TLibDecoder/TDecTop.h"
#include "TLibCommon/TComSlice.h"

using namespace std;

void Dec_Stat_Init()
{
    // Normalized Gaussian mask design, 11*11, s.d. 1.5
    int x, y;
    Double coeffSum = 0.0;
    for (y = 0; y < WEIGHTING_SIZE; y++)
    {
        for (x = 0; x < WEIGHTING_SIZE; x++)
        {
            decTopStat.weights[y][x] = exp(-((y - WEIGHTING_MID_TAP) * (y - WEIGHTING_MID_TAP) + (x - WEIGHTING_MID_TAP) * (x - WEIGHTING_MID_TAP)) / (WEIGHTING_MID_TAP - 0.5));
            coeffSum += decTopStat.weights[y][x];
        }
    }

    for (y = 0; y < WEIGHTING_SIZE; y++)
    {
        for (x = 0; x < WEIGHTING_SIZE; x++)
        {
            decTopStat.weights[y][x] /= coeffSum;
        }
    }
}

void Dec_Stat_DeInit()
{
    if(decTopStat.m_source_fp != nullptr)
        fclose(decTopStat.m_source_fp);
}

#define SSIM_EN     0
#define YUV_FACTOR  6
#if CALC_FLICK
vector<Pel> prvRecY;
vector<UChar> prvSrcY;
#endif
void Dec_Calc_PSNR_SSIM(TComPic *pcPic, int* decFrameLimit)
{
    int x, y, comp;
    TComPicYuv* pPicRec       = pcPic->getPicYuvRec();
    const int confLeft        = 0; // pcPic->getConformanceWindow().getWindowLeftOffset();
    const int confRight       = 0; // pcPic->getConformanceWindow().getWindowRightOffset();
    const int confTop         = 0; // pcPic->getConformanceWindow().getWindowTopOffset();
    const int confBottom      = 0; // pcPic->getConformanceWindow().getWindowBottomOffset();
    
    int widthY  = pPicRec->getWidth(COMPONENT_Y);
    int heightY = pPicRec->getHeight(COMPONENT_Y);
    const int strideY = pPicRec->getStride(COMPONENT_Y);
    const int planeOffsetY = (confLeft >> 0) + (confTop >> 0) * strideY;
    const int sizeY = widthY * heightY;
    const int sizeRec = strideY * heightY;

    int widthUV = pPicRec->getWidth(COMPONENT_Cb);
    int heightUV = pPicRec->getHeight(COMPONENT_Cb);
    const int strideUV = pPicRec->getStride(COMPONENT_Cb);
    const int planeOffsetUV = (confLeft >> 1) + (confTop >> 1) * strideUV;
    const int sizeUV = widthUV * heightUV;

    DecFrameStat framestat;
    Pel *imgRecY, *imgRecUV[2];

    vector<UChar> imgSrcY(sizeY, 0);
    vector<UChar> imgSrcUV[2];

    // Calc Frame Avg QP
    const int ctu_num = pcPic->getNumberOfCtusInFrame();
    int QPSum = 0;
    int QPCnt = 0;
    framestat.frame_maxQP = 0;
    framestat.frame_minQP = 0xff;
    for (x = 0; x < ctu_num; x++)
    {
        int qp_val = decTopStat.CTUStat[x].m_avgQP;
        if (qp_val != 0)
        {
            QPCnt++;
            QPSum += qp_val;
            framestat.frame_maxQP = max(decTopStat.CTUStat[x].m_maxQP, framestat.frame_maxQP);
            framestat.frame_minQP = min(decTopStat.CTUStat[x].m_minQP, framestat.frame_minQP);
        }
    }
    int QPAvg = (QPCnt == 0) ? pcPic->getSlice(0)->getSliceQp() : ((QPSum + (QPCnt >> 1)) / QPCnt);


    imgRecY     = pPicRec->getAddr(COMPONENT_Y) + planeOffsetY;
    imgRecUV[0] = pPicRec->getAddr(COMPONENT_Cb) + planeOffsetUV;
    imgRecUV[1] = pPicRec->getAddr(COMPONENT_Cr) + planeOffsetUV;

    imgSrcUV[0].resize(sizeUV, UChar(0));
    imgSrcUV[1].resize(sizeUV, UChar(0));

#if CALC_FLICK
    Int64 cntFlicker = 0;
    if (prvRecY.size() == 0)
        prvRecY.resize(sizeRec);
    if (prvSrcY.size() == 0)
        prvSrcY.resize(sizeY);
#endif

    framestat.frame_POC = pcPic->getPOC();
    framestat.frame_SliceQP = pcPic->getSlice(0)->getSliceQp();
    framestat.frame_AvgQP = QPAvg;
    framestat.slice_type = pcPic->getSlice(0)->getSliceType();
    framestat.frameBits = decTopStat.m_frameBits;
    framestat.SAOBits = decTopStat.m_SAOBits;
    framestat.HdrBits = decTopStat.m_HdrBits;
    framestat.DataBits = decTopStat.m_DBits;
    framestat.Flicker = 0;
    decTopStat.m_frameBits = 0;
    decTopStat.m_SAOBits = 0;
    decTopStat.m_HdrBits = 0;
    decTopStat.m_DBits = 0;

    // for YUV 420 3-plane
    if(decTopStat.m_source_fp != nullptr)
    {
        int chk_size;
        chk_size  = fread(&imgSrcY[0], sizeof(UChar), sizeY, decTopStat.m_source_fp);
        chk_size += fread(&imgSrcUV[0][0], sizeof(UChar), sizeUV, decTopStat.m_source_fp);
        chk_size += fread(&imgSrcUV[1][0], sizeof(UChar), sizeUV, decTopStat.m_source_fp);
        if (chk_size != (sizeY + 2 * sizeUV)*sizeof(UChar))
        {
        #if LOOP_INPUT
            _fseeki64(decTopStat.m_source_fp, 0, SEEK_SET);
            chk_size = fread(&imgSrcY[0], sizeof(UChar), sizeY, decTopStat.m_source_fp);
            chk_size += fread(&imgSrcUV[0][0], sizeof(UChar), sizeUV, decTopStat.m_source_fp);
            chk_size += fread(&imgSrcUV[1][0], sizeof(UChar), sizeUV, decTopStat.m_source_fp);
        #else
            memset(&imgSrcY[0], 0, sizeof(UChar) * imgSrcY.size());
            memset(&imgSrcUV[0][0], 0, sizeof(UChar) * imgSrcUV[0].size());
            memset(&imgSrcUV[1][0], 0, sizeof(UChar) * imgSrcUV[1].size());
            // discard this frame and frames hereafter
            decTopStat.m_decFrameNum--;
            *decFrameLimit = decTopStat.m_decFrameNum;
        #endif
        }
    }
    else
    {
        for (comp = 0; comp < 4; comp++)
        {
            framestat.m_decPSNR[comp] = 0;
            framestat.m_decSSIM[comp] = 0;
            //framestat.m_decMSE[comp] = 0;
        }
        decTopStat.FrameStat.push_back(framestat);
        return;
    }

    Int64 tmp = 0;
    float diff_avg;
    float diff_yuv;
    float psnr;
    const float psnr_max = 60.0;
    const float eps = (float)1e-10;
    const int maxValue = 255;
    widthY  -= (confLeft + confRight);  // conf_window should be zero
    heightY -= (confTop + confBottom);

    const TComSPS sps = pcPic->getPicSym()->getSPS();
    const int ctu_size = sps.getMaxCUWidth();
    const int ctu_sizeUV = ctu_size / 2;
    const int frame_width_ctu = (widthY + ctu_size - 1) / ctu_size;
    //const int frame_height_ctu = (heightY + ctu_size - 1) / ctu_size;
    //const int ctu_num = frame_width_ctu * frame_height_ctu;
    Int64 *ctuDiff = new Int64[ctu_num];
    memset(ctuDiff, 0, sizeof(Int64)* ctu_num);
    
    // Y
    for(y=0; y<heightY; y++)
    {
        for (x = 0; x < widthY; x++)
        {
            const int ctuX = x / ctu_size;
            const int ctuY = y / ctu_size;
            const int ctuIdx = ctuY * frame_width_ctu + ctuX;
            const int recV = imgRecY[y * strideY + x];
            const int srcV = imgSrcY[y * widthY + x];
            const int diff = recV - srcV;
            const Int64 diff2 = (Int64)diff * diff;
            tmp += diff2;
            ctuDiff[ctuIdx] += diff2;

            //fwrite(&imgSrcY[y * widthY + x], 1, 1, fpS);
            //fwrite(&imgRecY[y * strideY + x], 1, 1, fpR);
        #if CALC_FLICK
            if (framestat.slice_type == 2 && decTopStat.m_decFrameNum > 0)
            {
                const int PrecV = prvRecY[y * strideY + x];
                const int PsrcV = prvSrcY[y * widthY + x];
                cntFlicker += max(0, abs(PrecV - recV) - abs(PsrcV - srcV));
            }
        #endif
        }
    }

    diff_avg = (float)((double)tmp / (widthY * heightY));
    psnr = min(10 * log10(maxValue * maxValue / max(diff_avg, eps)), psnr_max);
    //psnr = (tmp == 0) ? psnr_max : (log10((double)((maxValue * maxValue * widthY * heightY) / tmp)) * 10.0);
    diff_yuv = psnr * YUV_FACTOR;
    framestat.m_decPSNR[0] = psnr;
    framestat.m_decMSE[0] = diff_avg;
#if CALC_FLICK
    if (framestat.slice_type == 2 && decTopStat.m_decFrameNum > 0)
        framestat.Flicker = cntFlicker / (widthY / 32 * heightY / 32);
    if (decTopStat.m_source_fp != nullptr)
    {
        memcpy(&prvRecY[0], &imgRecY[0], sizeof(Pel)* sizeRec);
        memcpy(&prvSrcY[0], &imgSrcY[0], sizeof(UChar) * sizeY);
    }
#endif

#if CTU_STAT_EN
    // ctu PSNR-Y
    for (x = 0; x < ctu_num; x++)
    {
        diff_avg = (float)((double)ctuDiff[x] / (ctu_size * ctu_size));
        psnr = min(10 * log10(maxValue * maxValue / max(diff_avg, eps)), psnr_max);
        decTopStat.CTUStat[x].m_ctuPSNR[0] = psnr;
        decTopStat.CTUStat[x].m_ctuMSE[0] = diff_avg;
        ctuDiff[x] = 0;
    }
#endif

    widthUV  -= ((confLeft>>1) + (confRight>>1));  // conf_window should be zero
    heightUV -= ((confTop>>1) + (confBottom>>1));
    const int frame_width_ctuUV = (widthUV + ctu_sizeUV - 1) / ctu_sizeUV;
    // UV
    for (comp = 0; comp < 2; comp++)
    {
        tmp = 0;
        for (y = 0; y < heightUV; y++)
        {
            for (x = 0; x < widthUV; x++)
            {
                const int ctuX = x / ctu_sizeUV;
                const int ctuY = y / ctu_sizeUV;
                const int ctuIdx = ctuY * frame_width_ctuUV + ctuX;
                const int recV = imgRecUV[comp][y * strideUV + x];
                const int srcV = imgSrcUV[comp][y * widthUV + x];
                const int diff = recV - srcV;
                const Int64 diff2 = (Int64)diff * diff;
                tmp += diff2;
                ctuDiff[ctuIdx] += diff2;
            }
        }
        diff_avg = (float)((double)tmp / (widthUV * heightUV));
        psnr = min(10 * log10(maxValue * maxValue / max(diff_avg, eps)), psnr_max);
        diff_yuv += psnr;
        framestat.m_decPSNR[1 + comp] = psnr;
        framestat.m_decMSE[1 + comp] = diff_avg;
#if CTU_STAT_EN
        // ctu PSNR-UV
        for (x = 0; x < ctu_num; x++)
        {
            diff_avg = (float)((double)ctuDiff[x] / (ctu_sizeUV * ctu_sizeUV));
            psnr = min(10 * log10(maxValue * maxValue / max(diff_avg, eps)), psnr_max);
            decTopStat.CTUStat[x].m_ctuPSNR[1 + comp] = psnr;
            decTopStat.CTUStat[x].m_ctuMSE[1 + comp] = diff_avg;
            ctuDiff[x] = 0;
        }
#endif
    }

    // PSNR-YUV
    //diff_yuv /= (YUV_FACTOR+2);
    //psnr = min(10 * log10(maxValue * maxValue / max(diff_yuv, eps)), psnr_max);
    psnr = diff_yuv / (YUV_FACTOR + 2);
    framestat.m_decPSNR[3] = psnr;

#if SSIM_EN
    // SSIM, this is HM way, while vmaf use downsampled reference
    {
        const Double c1 = (0.01 * maxValue) * (0.01 * maxValue);
        const Double c2 = (0.03 * maxValue) * (0.03 * maxValue);

        int blocksPerRow = widthY - WEIGHTING_SIZE + 1;
        int blocksPerColumn = heightY - WEIGHTING_SIZE + 1;
        int totalBlocks = blocksPerRow * blocksPerColumn;

        Double meanSSIM = 0.0;
        Double meanSSIM_yuv = 0.0;

        // SSIM-Y
        for (int blockIndexY = 0; blockIndexY < blocksPerColumn; blockIndexY++)
        {
            for (int blockIndexX = 0; blockIndexX < blocksPerRow; blockIndexX++)
            {
                Double muOrg = 0.0;
                Double muRec = 0.0;
                Double muOrigSqr = 0.0;
                Double muRecSqr = 0.0;
                Double muOrigMultRec = 0.0;

                for (y = 0; y < WEIGHTING_SIZE; y++)
                {
                    for (x = 0; x < WEIGHTING_SIZE; x++)
                    {
                        const Double gaussianWeight = decTopStat.weights[y][x];
                        //const Double recV = imgRecY[(blockIndexY + y) * strideY + (blockIndexX + x)];
                        //const Double srcV = imgSrcY[(blockIndexY + y) * widthY + (blockIndexX + x)];
                        const Pel recV = imgRecY[(blockIndexY + y) * strideY + (blockIndexX + x)];
                        const Pel srcV = imgSrcY[(blockIndexY + y) * widthY + (blockIndexX + x)];

                        Double recV_d = recV * gaussianWeight;
                        Double srcV_d = srcV * gaussianWeight;

                        muOrg += srcV_d;
                        muRec += recV_d;
                        muOrigSqr += srcV * srcV_d;
                        muRecSqr += recV * recV_d;
                        muOrigMultRec += srcV * recV_d;
                    }
                }

                const Double muOrg2 = muOrg * muOrg;
                const Double muRec2 = muRec * muRec;
                const Double muOrgRec = muOrg * muRec;
                const Double sigmaSqrOrig = muOrigSqr - muOrg2;
                const Double sigmaSqrRec = muRecSqr - muRec2;
                const Double sigmaOrigRec = muOrigMultRec - muOrgRec;

                Double blockSSIMVal = ((2.0 * sigmaOrigRec + c2) / (sigmaSqrOrig + sigmaSqrRec + c2));
                blockSSIMVal *= (2.0 * muOrgRec + c1) / (muOrg2 + muRec2 + c1);

                meanSSIM += blockSSIMVal;
            }
        }
        meanSSIM /= totalBlocks;
        meanSSIM_yuv += meanSSIM * YUV_FACTOR;
        framestat.m_decSSIM[0] = (float)meanSSIM;

        // SSIM-UV
        blocksPerRow = widthUV - WEIGHTING_SIZE + 1;
        blocksPerColumn = heightUV - WEIGHTING_SIZE + 1;
        totalBlocks = blocksPerRow * blocksPerColumn;
        for (comp = 0; comp < 2; comp++)
        {
            meanSSIM = 0.0;
            for (int blockIndexY = 0; blockIndexY < blocksPerColumn; blockIndexY++)
            {
                for (int blockIndexX = 0; blockIndexX < blocksPerRow; blockIndexX++)
                {
                    Double muOrg = 0.0;
                    Double muRec = 0.0;
                    Double muOrigSqr = 0.0;
                    Double muRecSqr = 0.0;
                    Double muOrigMultRec = 0.0;

                    for (y = 0; y < WEIGHTING_SIZE; y++)
                    {
                        for (x = 0; x < WEIGHTING_SIZE; x++)
                        {
                            const Double gaussianWeight = decTopStat.weights[y][x];
                            //const Double recV = imgRecUV[comp][(blockIndexY + y) * strideUV + (blockIndexX + x)];
                            //const Double srcV = imgSrcUV[comp][(blockIndexY + y) * widthUV + (blockIndexX + x)];
                            const Pel recV = imgRecUV[comp][(blockIndexY + y) * strideUV + (blockIndexX + x)];
                            const Pel srcV = imgSrcUV[comp][(blockIndexY + y) * widthUV + (blockIndexX + x)];

                            Double recV_d = recV * gaussianWeight;
                            Double srcV_d = srcV * gaussianWeight;

                            muOrg += srcV_d;
                            muRec += recV_d;
                            muOrigSqr += srcV * srcV_d;
                            muRecSqr += recV * recV_d;
                            muOrigMultRec += srcV * recV_d;
                        }
                    }

                    const Double muOrg2 = muOrg * muOrg;
                    const Double muRec2 = muRec * muRec;
                    const Double muOrgRec = muOrg * muRec;
                    const Double sigmaSqrOrig = muOrigSqr - muOrg2;
                    const Double sigmaSqrRec = muRecSqr - muRec2;
                    const Double sigmaOrigRec = muOrigMultRec - muOrgRec;

                    Double blockSSIMVal = ((2.0 * sigmaOrigRec + c2) / (sigmaSqrOrig + sigmaSqrRec + c2));
                    blockSSIMVal *= (2.0 * muOrgRec + c1) / (muOrg2 + muRec2 + c1);

                    meanSSIM += blockSSIMVal;
                }
            }
            meanSSIM /= totalBlocks;
            meanSSIM_yuv += meanSSIM;
            framestat.m_decSSIM[1 + comp] = (float)meanSSIM;
        }
        meanSSIM_yuv /= (YUV_FACTOR+2);
        framestat.m_decSSIM[3] = (float)meanSSIM_yuv;
    }
#else
    framestat.m_decSSIM[0] = 0;
    framestat.m_decSSIM[1] = 0;
    framestat.m_decSSIM[2] = 0;
    framestat.m_decSSIM[3] = 0;
#endif
    decTopStat.FrameStat.push_back(framestat);

    delete[] ctuDiff;
}


#define ONLY_YUV    1
#define FRM_STAT    1
#define toK(x)  ((float)(x)/1000)
void Print_RC_Stat()
{
    FILE* fp_rc = fopen("rc_stat.csv", "w");
    int i, nb;
#if ONLY_YUV
    const int comp = 3;
#else   
    int comp = 0;
#endif

    float avg_psnr[4] = { 0, 0, 0, 0 };
    float max_psnr[4] = { 0, 0, 0, 0 };
    float min_psnr[4] = { 255, 255, 255, 255 };
    double var_psnr[4] = { 0, 0, 0, 0 };
    double psnr2[4] = { 0, 0, 0, 0 };
    
    float avg_ssim[4] = { 0, 0, 0, 0 };
    float max_ssim[4] = { 0, 0, 0, 0 };
    float min_ssim[4] = { 255, 255, 255, 255 };
    double var_ssim[4] = { 0, 0, 0, 0 };
    double ssim2[4] = { 0, 0, 0, 0 };
#if CALC_FLICK
    int avgFlick = 0;
    int cntFlick = 0;
    decTopStat.m_avgFlick = 0;
#endif

    if (fp_rc == nullptr)
    {
        printf("RC: open stat file failed!\n");
        return;
    }

    // Find Average/Max/Min PSNR/SSIM
    const int frame_num = decTopStat.m_decFrameNum;
    for (nb = 0; nb < frame_num; nb++)
    {
#if !ONLY_YUV
        for (comp = 0; comp < 4; comp++)
#endif
        {
            avg_psnr[comp] += decTopStat.FrameStat[nb].m_decPSNR[comp];
            psnr2[comp] += Double(decTopStat.FrameStat[nb].m_decPSNR[comp]) * decTopStat.FrameStat[nb].m_decPSNR[comp];
            if (max_psnr[comp] < decTopStat.FrameStat[nb].m_decPSNR[comp])
                max_psnr[comp] = decTopStat.FrameStat[nb].m_decPSNR[comp];
            if (min_psnr[comp] > decTopStat.FrameStat[nb].m_decPSNR[comp])
                min_psnr[comp] = decTopStat.FrameStat[nb].m_decPSNR[comp];

            avg_ssim[comp] += decTopStat.FrameStat[nb].m_decSSIM[comp];
            ssim2[comp] += Double(decTopStat.FrameStat[nb].m_decSSIM[comp]) * decTopStat.FrameStat[nb].m_decSSIM[comp];
            if (max_ssim[comp] < decTopStat.FrameStat[nb].m_decSSIM[comp])
                max_ssim[comp] = decTopStat.FrameStat[nb].m_decSSIM[comp];
            if (min_ssim[comp] > decTopStat.FrameStat[nb].m_decSSIM[comp])
                min_ssim[comp] = decTopStat.FrameStat[nb].m_decSSIM[comp];
        }
        decTopStat.m_maxQP = max(decTopStat.FrameStat[nb].frame_maxQP, decTopStat.m_maxQP);
        decTopStat.m_minQP = min(decTopStat.FrameStat[nb].frame_minQP, decTopStat.m_minQP);
    #if CALC_FLICK
        if (decTopStat.FrameStat[nb].slice_type == 2 && nb > 0)
        {
            avgFlick += decTopStat.FrameStat[nb].Flicker;
            cntFlick++;
        }
    #endif
    }
    
    // Find Mean/Variance PSNR/SSIM
#if !ONLY_YUV
    for (comp = 0; comp < 4; comp++)
#endif
    {
        avg_psnr[comp] /= frame_num;
        avg_ssim[comp] /= frame_num;

        var_psnr[comp] = (psnr2[comp] / frame_num) - Double(avg_psnr[comp]) * avg_psnr[comp];
        var_psnr[comp] = sqrt(var_psnr[comp]);
        var_ssim[comp] = (ssim2[comp] / frame_num) - Double(avg_ssim[comp]) * avg_ssim[comp];
        var_ssim[comp] = sqrt(var_ssim[comp]);
    #if CALC_FLICK
        if (cntFlick != 0)
            decTopStat.m_avgFlick = avgFlick / cntFlick;
    #endif
    }

    // Find Average Rate in GOPs
    int gop_num = max(min((int)decTopStat.m_GOPnum, frame_num), 1);
    const int GOPs = max(frame_num / gop_num, 1);
    const int target_rate = decTopStat.m_targetRate;
    vector<int> GOP_RATE(GOPs+1, 0);
    vector<int> GOP_MEAN(GOPs+1, 0);
    
    int avg_frm_rate = 0;
    int max_frm_rate = 0;
    int min_frm_rate = MAX_INT;
    Int64 var_frm_rate = 0;
    Int64 var_frm_rate2 = 0;
    int avg_gop_rate = 0;
    int max_gop_rate = 0;
    int min_gop_rate = MAX_INT;
    Int64 var_gop_rate = 0;
    Int64 var_gop_rate2 = 0;
    Int64 all_gop_rate = 0;

    for (nb = 0; nb < GOPs; nb++)
    {
        int gop_rate = 0;
        for (i = 0; i < gop_num; i++)
        {
            const int framebits = decTopStat.FrameStat[nb * gop_num + i].frameBits;
            gop_rate += framebits;
            var_frm_rate2 += Int64(framebits) * framebits;
#if FRM_STAT
            if (max_frm_rate < framebits)
                max_frm_rate = framebits;
            if (min_frm_rate > framebits)
                min_frm_rate = framebits;
#endif
        }
        GOP_RATE[nb] = gop_rate;
        gop_rate = (gop_rate + gop_num - 1) / gop_num;
        GOP_MEAN[nb] = gop_rate;

        all_gop_rate += GOP_RATE[nb];
        var_gop_rate2 += Int64(GOP_RATE[nb]) * GOP_RATE[nb];
        if (max_gop_rate < GOP_RATE[nb])
            max_gop_rate = GOP_RATE[nb];
        if (min_gop_rate > GOP_RATE[nb])
            min_gop_rate = GOP_RATE[nb];
    }
    avg_gop_rate = (int)((all_gop_rate + GOPs - 1) / GOPs);
    avg_frm_rate = (avg_gop_rate + gop_num - 1) / gop_num;

    var_gop_rate = ((var_gop_rate2 + GOPs - 1) / GOPs) - Int64(avg_gop_rate) * avg_gop_rate;
    var_gop_rate = (Int64)(sqrt((long double)var_gop_rate) + 0.5);

#if FRM_STAT
    var_frm_rate = (var_frm_rate2 + GOPs * gop_num - 1) / (GOPs * gop_num) - Int64(avg_frm_rate) * avg_frm_rate;
    var_frm_rate = (Int64)(sqrt((long double)var_frm_rate) + 0.5);
    /*for (nb = 0; nb < GOPs; nb++)
    {
        Int64 tmpv = 0;
        for (i = 0; i < gop_num; i++)
        {
            tmpi = decTopStat.FrameStat[nb * gop_num + i].frameBits - avg_frm_rate;
            tmpv += tmpi * tmpi;
        }
        var_frm_rate += (tmpv + GOPs * gop_num - 1) / (GOPs * gop_num);
    }
    var_frm_rate = (Int64)(sqrt(var_frm_rate) + 0.5);*/
#endif

    // Find CPB fullness
    gop_num = decTopStat.m_GOPnum;
    const int BufferingRate = (target_rate + gop_num - 1) / gop_num;
    int CPBLoBound = BufferingRate;
    int CPBUpBound;
    int CPBRecover;
    int CPBState;
    int gop_bits = target_rate;
    vector<int> FRM_CPB(frame_num, 0);
    vector<int> CPBover(frame_num, 0);
    vector<int> CPBuder(frame_num, 0);
    int CPBover_sum = 0;
    int CPBuder_sum = 0;
    float fullFactor = 0.9f;

    if (decTopStat.m_targetCpbSize == 0)
    {
        // CR = 0.8/1.5, for level 5/6, MinCR=6 or 8, there will be 3 to 4 bitstream unit buffers
        //int defaultCPB = (int)(decTopStat.m_frameWidth * decTopStat.m_frameHeight * 0.8 + 0.5);
        int defaultCPB = BufferingRate * 25;
        decTopStat.m_targetCpbSize = min(defaultCPB, decTopStat.m_MaxCpbSize);
    }
    else if (decTopStat.m_targetCpbSize == 1)
    {
        decTopStat.m_targetCpbSize = decTopStat.m_MaxCpbSize;
        fullFactor = 0.9f;
    }
    else if (decTopStat.m_targetCpbSize > decTopStat.m_MaxCpbSize)
    {
        decTopStat.m_targetCpbSize = decTopStat.m_MaxCpbSize;
    }
    CPBUpBound = (int)(decTopStat.m_targetCpbSize * 0.9 + 0.5);
    CPBRecover = (int)(decTopStat.m_targetCpbSize * 0.5 + 0.5);

    nb = 0;
    CPBState = (int)(decTopStat.m_targetCpbSize * fullFactor + 0.5);
    CPBState += (BufferingRate - decTopStat.FrameStat[nb].frameBits);
    gop_bits -= decTopStat.FrameStat[nb].frameBits;
    if (CPBState > CPBUpBound)
    {
        CPBover[nb] = 1;
        CPBover_sum++;
        // overflow, wait CPB to drain untill 0.5*MAX_CPB_BUF
        CPBState = CPBRecover;
    }
    if ((CPBState - BufferingRate) < CPBLoBound)
    {
        CPBuder[nb] = -1;
        CPBuder_sum++;
        // underflow, wait CPB to fill untill 0.5*MAX_CPB_BUF
        CPBState = CPBRecover;
    }
    FRM_CPB[nb] = CPBState;

    for (nb = 1; nb < frame_num; nb++)
    {
        int next_frm_pos = (nb + 1) % gop_num;
        int left_frm = gop_num - next_frm_pos;
        gop_bits -= decTopStat.FrameStat[nb].frameBits;
        if (next_frm_pos == 0)
        {
            gop_bits += target_rate;
        }
        CPBLoBound = (gop_bits + left_frm - 1) / left_frm;
        if (CPBLoBound < 100)
        {
            CPBLoBound = 100;   // at least allocate 100 bits for one picture
        }
        CPBLoBound = int(CPBLoBound * 0.1 + BufferingRate * 0.9);

        CPBState += (BufferingRate - decTopStat.FrameStat[nb].frameBits);
        if (CPBState > CPBUpBound)
        {
            CPBover[nb] = 1;
            CPBover_sum++;
            // overflow, wait CPB to drain untill 0.5*MAX_CPB_BUF
            CPBState = CPBRecover;
        }
        if ((CPBState - BufferingRate) < CPBLoBound)
        {
            CPBuder[nb] = -1;
            CPBuder_sum++;
            // underflow, wait CPB to fill untill 0.5*MAX_CPB_BUF
            CPBState = CPBRecover;
        }
        FRM_CPB[nb] = CPBState;
    }

    // AVBR file size check
    int estFileSize = GOPs * target_rate;
    float file_size_diff = (float)(all_gop_rate - estFileSize) * 100 / estFileSize;

#if RCSTAT_COM == 0
    // Print Frame Info
    fprintf(fp_rc, "Bitstream: %s\n", decTopStat.streamName.c_str());
    fprintf(fp_rc, "Width x Height: %d x %d\n", decTopStat.m_frameWidth, decTopStat.m_frameHeight);
    fprintf(fp_rc, "Total Frame Counts: %d\n", frame_num);
    fprintf(fp_rc, "Total GOPs: %d\n", GOPs);
    fprintf(fp_rc, "CPB Max Size: %.2f Kbits\n", toK(decTopStat.m_targetCpbSize));
    fprintf(fp_rc, "Target Rate of GOP(%d): %.2f Kbits\n", gop_num, toK(target_rate));
    fprintf(fp_rc, "Total CPB Overflow Counts = %d\n", CPBover_sum);
    fprintf(fp_rc, "Total CPB Underflow Counts = %d\n", CPBuder_sum);
    fprintf(fp_rc, "AVBR Estimated File Size: %.2f Kbits\n", toK(estFileSize));
    fprintf(fp_rc, "AVBR Actual File Size: %.2f Kbits\n", toK(all_gop_rate));
    fprintf(fp_rc, "AVBR File Size Diff: %+2.2f%%\n", file_size_diff);
    fprintf(fp_rc, "\n");

    // Print Frame Average
    fprintf(fp_rc, ", MAX, MIN, AVG, STD\n");
    fprintf(fp_rc, "QP, %d, %d, 0, 0\n", decTopStat.m_maxQP, decTopStat.m_minQP);
#if !ONLY_YUV
    fprintf(fp_rc, "PSNR-Y, %2.2f, %2.2f, %2.2f, %2.2f\n", max_psnr[0], min_psnr[0], avg_psnr[0], (float)var_psnr[0]);
    fprintf(fp_rc, "PSNR-U, %2.2f, %2.2f, %2.2f, %2.2f\n", max_psnr[1], min_psnr[1], avg_psnr[1], (float)var_psnr[1]);
    fprintf(fp_rc, "PSNR-V, %2.2f, %2.2f, %2.2f, %2.2f\n", max_psnr[2], min_psnr[2], avg_psnr[2], (float)var_psnr[2]);
#endif
    fprintf(fp_rc, "PSNR-YUV, %2.2f, %2.2f, %2.2f, %2.2f, dB\n", max_psnr[3], min_psnr[3], avg_psnr[3], (float)var_psnr[3]);
#if !ONLY_YUV
    fprintf(fp_rc, "SSIM-Y, %1.4f, %1.4f, %1.4f, %1.4f\n", max_ssim[0], min_ssim[0], avg_ssim[0], (float)var_ssim[0]);
    fprintf(fp_rc, "SSIM-U, %1.4f, %1.4f, %1.4f, %1.4f\n", max_ssim[1], min_ssim[1], avg_ssim[1], (float)var_ssim[1]);
    fprintf(fp_rc, "SSIM-V, %1.4f, %1.4f, %1.4f, %1.4f\n", max_ssim[2], min_ssim[2], avg_ssim[2], (float)var_ssim[2]);
#endif
    fprintf(fp_rc, "SSIM-YUV, %1.4f, %1.4f, %1.4f, %1.4f\n", max_ssim[3], min_ssim[3], avg_ssim[3], (float)var_ssim[3]);
    fprintf(fp_rc, "GOP Rate, %.2f, %.2f, %.2f, %.2f, KBits\n", toK(max_gop_rate), toK(min_gop_rate), toK(avg_gop_rate), toK(var_gop_rate));
#if FRM_STAT
    fprintf(fp_rc, "FRM Rate, %.2f, %.2f, %.2f, %.2f, KBits\n\n", toK(max_frm_rate), toK(min_frm_rate), toK(avg_frm_rate), toK(var_frm_rate));
#endif
#endif

    // Print Frame Details
    nb = 0;
#if CHK_TX_QP
    fprintf(fp_rc, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n",
            "FrameNo.", "POC", "Slice", "SliceQP", "AvgQP", "MaxQP", "MinQP",
            "PSNR-YUV", "PSNR-Prev", "PSNR-Avg",
            "SSIM-YUV", "SSIM-Prev", "SSIM-Avg",
            "Bits", "Bits-Avg", "Bits-FRMAvg", "GOP-GOPAvg", /*GOP-Target, */"CPB Full", "CPB CHK", "SAOBits", "HdrBits", "DataBits"
    );
    fprintf(fp_rc, "%d, %d, %s, %d, %d, %d, %d, %2.2f, %+2.2f, %+2.2f, %1.4f, %+1.4f, %+1.4f, %d, %+d, %+d, ",
            nb,
            decTopStat.FrameStat[nb].frame_POC,
            (decTopStat.FrameStat[nb].slice_type == 2) ? "I" : (decTopStat.FrameStat[nb].slice_type == 1) ? "P" : "B",
            decTopStat.FrameStat[nb].frame_SliceQP,
            decTopStat.FrameStat[nb].frame_AvgQP,
            decTopStat.FrameStat[nb].frame_maxQP,
            decTopStat.FrameStat[nb].frame_minQP,
            decTopStat.FrameStat[nb].m_decPSNR[3], 0.0, decTopStat.FrameStat[nb].m_decPSNR[3] - avg_psnr[3],
            decTopStat.FrameStat[nb].m_decSSIM[3], 0.0, decTopStat.FrameStat[nb].m_decSSIM[3] - avg_ssim[3],
            /*toK*/(decTopStat.FrameStat[nb].frameBits), /*toK*/(decTopStat.FrameStat[nb].frameBits - GOP_MEAN[nb / gop_num]), /*toK*/(decTopStat.FrameStat[nb].frameBits - avg_frm_rate)
    );
#else
    fprintf(fp_rc, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n",
            "FrameNo.", "POC", "Slice", "SliceQP", "AvgQP",
            "PSNR-YUV", "PSNR-Prev", "PSNR-Avg",
            "SSIM-YUV", "SSIM-Prev", "SSIM-Avg",
            "Bits", "Bits-Avg", "Bits-FRMAvg", "GOP-GOPAvg", /*GOP-Target, */"CPB Full", "CPB CHK", "SAOBits", "HdrBits", "DataBits"
    );
    fprintf(fp_rc, "%d, %d, %s, %d, %d, %2.2f, %+2.2f, %+2.2f, %1.4f, %+1.4f, %+1.4f, %d, %+d, %+d, ",
            nb,
            decTopStat.FrameStat[nb].frame_POC,
            (decTopStat.FrameStat[nb].slice_type == 2) ? "I" : (decTopStat.FrameStat[nb].slice_type == 1) ? "P" : "B",
            decTopStat.FrameStat[nb].frame_SliceQP,
            decTopStat.FrameStat[nb].frame_AvgQP,
            decTopStat.FrameStat[nb].m_decPSNR[3], 0.0, decTopStat.FrameStat[nb].m_decPSNR[3] - avg_psnr[3],
            decTopStat.FrameStat[nb].m_decSSIM[3], 0.0, decTopStat.FrameStat[nb].m_decSSIM[3] - avg_ssim[3],
            /*toK*/(decTopStat.FrameStat[nb].frameBits), /*toK*/(decTopStat.FrameStat[nb].frameBits - GOP_MEAN[nb / gop_num]), /*toK*/(decTopStat.FrameStat[nb].frameBits - avg_frm_rate)
    );
#endif
    if ((nb % gop_num) == 0)
        fprintf(fp_rc, "%+d", /*toK*/(GOP_RATE[nb / gop_num] - avg_gop_rate));
    fprintf(fp_rc, ",%3.1f%%, ", (float)FRM_CPB[nb]*100/decTopStat.m_targetCpbSize);
    if(CPBover[nb])
        fprintf(fp_rc, "%d", CPBover[nb]);
    else if(CPBuder[nb])
        fprintf(fp_rc, "%d", CPBuder[nb]);

    // append
    fprintf(fp_rc, ",%d, %d, %d", decTopStat.FrameStat[nb].SAOBits, decTopStat.FrameStat[nb].HdrBits, decTopStat.FrameStat[nb].DataBits);

    fprintf(fp_rc, "\n");
    printf("Frame %d PSNR-Y: %2.2f, PSNR-U: %2.2f, PSNR-V: %2.2f\n", nb, decTopStat.FrameStat[nb].m_decPSNR[0], decTopStat.FrameStat[nb].m_decPSNR[1], decTopStat.FrameStat[nb].m_decPSNR[2]);
    
    for (nb = 1; nb < frame_num; nb++)
    {
#if CHK_TX_QP
        fprintf(fp_rc, "%d, %d, %s, %d, %d, %d, %d, %2.2f, %+2.2f, %+2.2f, %1.4f, %+1.4f, %+1.4f, %d, %+d, %+d, ",
                nb,
                decTopStat.FrameStat[nb].frame_POC,
                (decTopStat.FrameStat[nb].slice_type == 2) ? "I" : (decTopStat.FrameStat[nb].slice_type == 1) ? "P" : "B",
                decTopStat.FrameStat[nb].frame_SliceQP,
                decTopStat.FrameStat[nb].frame_AvgQP,
                decTopStat.FrameStat[nb].frame_maxQP,
                decTopStat.FrameStat[nb].frame_minQP,
                decTopStat.FrameStat[nb].m_decPSNR[3], decTopStat.FrameStat[nb].m_decPSNR[3] - decTopStat.FrameStat[nb - 1].m_decPSNR[3], decTopStat.FrameStat[nb].m_decPSNR[3] - avg_psnr[3],
                decTopStat.FrameStat[nb].m_decSSIM[3], decTopStat.FrameStat[nb].m_decSSIM[3] - decTopStat.FrameStat[nb - 1].m_decSSIM[3], decTopStat.FrameStat[nb].m_decSSIM[3] - avg_ssim[3],
                /*toK*/(decTopStat.FrameStat[nb].frameBits), /*toK*/(decTopStat.FrameStat[nb].frameBits - GOP_MEAN[nb / gop_num]), /*toK*/(decTopStat.FrameStat[nb].frameBits - avg_frm_rate)
        );
#else
        fprintf(fp_rc, "%d, %d, %s, %d, %d, %2.2f, %+2.2f, %+2.2f, %1.4f, %+1.4f, %+1.4f, %d, %+d, %+d, ",
                nb,
                decTopStat.FrameStat[nb].frame_POC,
                (decTopStat.FrameStat[nb].slice_type == 2) ? "I" : (decTopStat.FrameStat[nb].slice_type == 1) ? "P" : "B",
                decTopStat.FrameStat[nb].frame_SliceQP,
                decTopStat.FrameStat[nb].frame_AvgQP,
                decTopStat.FrameStat[nb].m_decPSNR[3], decTopStat.FrameStat[nb].m_decPSNR[3] - decTopStat.FrameStat[nb - 1].m_decPSNR[3], decTopStat.FrameStat[nb].m_decPSNR[3] - avg_psnr[3],
                decTopStat.FrameStat[nb].m_decSSIM[3], decTopStat.FrameStat[nb].m_decSSIM[3] - decTopStat.FrameStat[nb - 1].m_decSSIM[3], decTopStat.FrameStat[nb].m_decSSIM[3] - avg_ssim[3],
                /*toK*/(decTopStat.FrameStat[nb].frameBits), /*toK*/(decTopStat.FrameStat[nb].frameBits - GOP_MEAN[nb / gop_num]), /*toK*/(decTopStat.FrameStat[nb].frameBits - avg_frm_rate)
        );
#endif
        if ((nb % gop_num) == 0)
            fprintf(fp_rc, "%+d", /*toK*/(GOP_RATE[nb / gop_num] - avg_gop_rate));
        fprintf(fp_rc, ",%3.1f%%, ", (float)FRM_CPB[nb] * 100 / decTopStat.m_targetCpbSize);
        if (CPBover[nb])
            fprintf(fp_rc, "%d", CPBover[nb]);
        else if (CPBuder[nb])
            fprintf(fp_rc, "%d", CPBuder[nb]);

        // append
        fprintf(fp_rc, ",%d, %d, %d", decTopStat.FrameStat[nb].SAOBits, decTopStat.FrameStat[nb].HdrBits, decTopStat.FrameStat[nb].DataBits);

        fprintf(fp_rc, "\n");

        printf("Frame %d PSNR-Y: %2.2f, PSNR-U: %2.2f, PSNR-V: %2.2f\n", nb, decTopStat.FrameStat[nb].m_decPSNR[0], decTopStat.FrameStat[nb].m_decPSNR[1], decTopStat.FrameStat[nb].m_decPSNR[2]);
    }

#if RCSTAT_COM == 1
    // Print Frame Info
    fprintf(fp_rc, "\n");
    fprintf(fp_rc, "Bitstream: %s\n", decTopStat.streamName.c_str());
    fprintf(fp_rc, "Width x Height: %d x %d\n", decTopStat.m_frameWidth, decTopStat.m_frameHeight);
    fprintf(fp_rc, "Total Frame Counts: %d\n", frame_num);
    fprintf(fp_rc, "Total GOPs: %d\n", GOPs);
    fprintf(fp_rc, "CPB Max Size: %.2f Kbits\n", toK(decTopStat.m_targetCpbSize));
    fprintf(fp_rc, "Target Rate of GOP(%d): %.2f Kbits\n", gop_num, toK(target_rate));
    fprintf(fp_rc, "Total CPB Overflow Counts = %d\n", CPBover_sum);
    fprintf(fp_rc, "Total CPB Underflow Counts = %d\n", CPBuder_sum);
    fprintf(fp_rc, "AVBR Estimated File Size: %.2f Kbits\n", toK(estFileSize));
    fprintf(fp_rc, "AVBR Actual File Size: %.2f Kbits\n", toK(all_gop_rate));
    fprintf(fp_rc, "AVBR File Size Diff: %+2.2f%%\n", file_size_diff);
    fprintf(fp_rc, "Average Flicker: %d\n", decTopStat.m_avgFlick);
    fprintf(fp_rc, "\n");

    // Print Frame Average
    fprintf(fp_rc, ", MAX, MIN, AVG, STD\n");
    fprintf(fp_rc, "QP, %d, %d, 0, 0\n", decTopStat.m_maxQP, decTopStat.m_minQP);
#if !ONLY_YUV
    fprintf(fp_rc, "PSNR-Y, %2.2f, %2.2f, %2.2f, %2.2f\n", max_psnr[0], min_psnr[0], avg_psnr[0], (float)var_psnr[0]);
    fprintf(fp_rc, "PSNR-U, %2.2f, %2.2f, %2.2f, %2.2f\n", max_psnr[1], min_psnr[1], avg_psnr[1], (float)var_psnr[1]);
    fprintf(fp_rc, "PSNR-V, %2.2f, %2.2f, %2.2f, %2.2f\n", max_psnr[2], min_psnr[2], avg_psnr[2], (float)var_psnr[2]);
#endif
    fprintf(fp_rc, "PSNR-YUV, %2.2f, %2.2f, %2.2f, %2.2f, dB\n", max_psnr[3], min_psnr[3], avg_psnr[3], (float)var_psnr[3]);
#if !ONLY_YUV
    fprintf(fp_rc, "SSIM-Y, %1.4f, %1.4f, %1.4f, %1.4f\n", max_ssim[0], min_ssim[0], avg_ssim[0], (float)var_ssim[0]);
    fprintf(fp_rc, "SSIM-U, %1.4f, %1.4f, %1.4f, %1.4f\n", max_ssim[1], min_ssim[1], avg_ssim[1], (float)var_ssim[1]);
    fprintf(fp_rc, "SSIM-V, %1.4f, %1.4f, %1.4f, %1.4f\n", max_ssim[2], min_ssim[2], avg_ssim[2], (float)var_ssim[2]);
#endif
    fprintf(fp_rc, "SSIM-YUV, %1.4f, %1.4f, %1.4f, %1.4f\n", max_ssim[3], min_ssim[3], avg_ssim[3], (float)var_ssim[3]);
    fprintf(fp_rc, "GOP Rate, %.2f, %.2f, %.2f, %.2f, KBits\n", toK(max_gop_rate), toK(min_gop_rate), toK(avg_gop_rate), toK(var_gop_rate));
#if FRM_STAT
    fprintf(fp_rc, "FRM Rate, %.2f, %.2f, %.2f, %.2f, KBits\n\n", toK(max_frm_rate), toK(min_frm_rate), toK(avg_frm_rate), toK(var_frm_rate));
#endif
#endif

#if CALC_FLICK
    prvRecY.resize(0);
    prvSrcY.resize(0);
#endif

    fclose(fp_rc);
}


void Dec_CTUStat_Init(const TComSPS* sps, const TComPPS* pps)
{
    if (decTopStat.m_isCTUInit)
        return;

    const int frame_width = sps->getPicWidthInLumaSamples();
    const int frame_height = sps->getPicHeightInLumaSamples();
    const int ctu_size = sps->getMaxCUWidth();

    const int frame_width_ctu = (frame_width + ctu_size - 1) / ctu_size;
    const int frame_height_ctu = (frame_height + ctu_size - 1) / ctu_size;
    const int ctu_num = frame_width_ctu * frame_height_ctu;

    if (sps == nullptr || pps == nullptr)
    {
        printf("RC_STAT: SPS/PPS is NULL\n");
        return;
    }

    system("mkdir ctu_info");

    decTopStat.CTUStat = new DecCTUStat[ctu_num];

    decTopStat.m_isCTUInit = 1;
}

void Dec_CTUStat_DeInit()
{
    if (decTopStat.m_isCTUInit)
    {
        delete[] decTopStat.CTUStat;
    }
}

//#define COMBINE64
void Print_CTU_Stat(TComPic* pcPic)
{
    FILE* fp;
    char file_name[30];
    
    const TComSPS sps = pcPic->getPicSym()->getSPS();
    //const TComPPS pps = pcPic->getPicSym()->getPPS();

    int x, y, k, sum_h, sum_d;
    const int frame_width = sps.getPicWidthInLumaSamples();
    const int frame_height = sps.getPicHeightInLumaSamples();
    const int ctu_size = sps.getMaxCUWidth();

    const int frame_width_ctu = (frame_width + ctu_size - 1) / ctu_size;
    const int frame_height_ctu = (frame_height + ctu_size - 1) / ctu_size;
    const int ctu_num = frame_width_ctu * frame_height_ctu;

    int ctu_part = 64 / ctu_size;
    int ctu_part2 = ctu_part * ctu_part;


    sprintf(file_name, "%s%d%s", "ctu_info/frame_", decTopStat.m_decFrameNum, ".csv");
    fp = fopen(file_name, "w");
        
    sum_h = 0;
    sum_d = 0;
    fprintf(fp, "CTU, POS-X, POS-Y, QP, H_Bit, D_Bits, PSNR-Y, PSNR-U, PSNR-V\n");

#ifdef COMBINE64
    for (y = 0; y < frame_height_ctu; y += ctu_part)
    {
        for (x = 0; x < frame_width_ctu; x += ctu_part)
        {
            int ctuNum = (y / ctu_part) * (frame_width_ctu / ctu_part) + (x / ctu_part);
            int ctuPosX = x * ctu_size;
            int ctuPosY = y * ctu_size;
            int AvgQP = 0;
            int HdrBits = 0;
            int DataBits = 0;
            float avgPSNR[3] = { 0, 0, 0 };
            for (k = 0; k < ctu_part2; k++)
            {
                int ctuIdx = (y + (k >> 1)) * frame_width_ctu + (x + (k & 0x1));
                AvgQP += decTopStat.CTUStat[ctuIdx].m_avgQP;
                HdrBits += decTopStat.CTUStat[ctuIdx].m_ctuHeaderBits;
                DataBits += decTopStat.CTUStat[ctuIdx].m_ctuDataBits;
                avgPSNR[0] += decTopStat.CTUStat[ctuIdx].m_ctuPSNR[0];
                avgPSNR[1] += decTopStat.CTUStat[ctuIdx].m_ctuPSNR[1];
                avgPSNR[2] += decTopStat.CTUStat[ctuIdx].m_ctuPSNR[2];

                sum_h += decTopStat.CTUStat[ctuIdx].m_ctuHeaderBits;
                sum_d += decTopStat.CTUStat[ctuIdx].m_ctuDataBits;
            }
            AvgQP /= ctu_part2;
            avgPSNR[0] /= ctu_part2;
            avgPSNR[1] /= ctu_part2;
            avgPSNR[2] /= ctu_part2;
            fprintf(fp, "%d, %03d, %03d, %d, %d, %d, %2.2f, %2.2f, %2.2f\n", ctuNum, ctuPosX, ctuPosY, AvgQP, HdrBits, DataBits, avgPSNR[0], avgPSNR[1], avgPSNR[2]);
        }
    }
#else
    for (x = 0; x < ctu_num; x++)
    {
        fprintf(fp, "%d, %03d, %03d, %d, %d, %d, %2.2f, %2.2f, %2.2f\n", x, (x % frame_width_ctu)*ctu_size, (x / frame_width_ctu)*ctu_size,
                decTopStat.CTUStat[x].m_avgQP,
                decTopStat.CTUStat[x].m_ctuHeaderBits,
                decTopStat.CTUStat[x].m_ctuDataBits,
                decTopStat.CTUStat[x].m_ctuPSNR[0],
                decTopStat.CTUStat[x].m_ctuPSNR[1],
                decTopStat.CTUStat[x].m_ctuPSNR[2]
        );
        sum_h += decTopStat.CTUStat[x].m_ctuHeaderBits;
        sum_d += decTopStat.CTUStat[x].m_ctuDataBits;
    }
#endif
    fprintf(fp, "Total H_Bits=%d\n", sum_h);
    fprintf(fp, "Total D_Bits=%d\n", sum_d);
    fprintf(fp, "Total Sum Bits=%d\n", sum_h + sum_d);
#if RC_STAT_EN
    fprintf(fp, "Total Frame Bits=%d\n", decTopStat.FrameStat[decTopStat.m_decFrameNum].frameBits);
#endif

    memset(&decTopStat.CTUStat[0], 0, sizeof(DecCTUStat) * ctu_num);

    fclose(fp);
}