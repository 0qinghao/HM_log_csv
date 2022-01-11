
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include "TLibEncoder/TEncTop.h"


void PreAnalyze_HAD(TComPic* curPic)
{
    TComPicYuv* pcPicYuv = curPic->getPicYuvOrg();
    Pel* curY = pcPicYuv->getAddr(COMPONENT_Y);
    const int width = pcPicYuv->getWidth(COMPONENT_Y);
    const int height = pcPicYuv->getHeight(COMPONENT_Y);
    const int stride = pcPicYuv->getStride(COMPONENT_Y);
    const int ctb_size = curPic->getSlice(0)->getSPS()->getMaxCUWidth();
    const int blk_num = ctb_size >> 3;
    const int frameno = encTopStat.m_encFrameNum;
    const int frame_sw = frameno & 0x1;

    int x, y, blkx, blky;
    int ctu_idx = 0;
    int tmp = 0;

    for (y = 0; y < height; y += ctb_size)
    {
        for (x = 0; x < width; x += ctb_size)
        {
            int blk8_idx = 0;
            int blk8_avg = 0;
            for (blky = 0; blky < ctb_size; blky+=8)  // frame size must align to 8
            {
                int y_pos = y + blky;
                if (y_pos >= height)
                    break;

                for (blkx = 0; blkx < ctb_size; blkx+=8)
                {
                    int k, i, j, jj;
                    int diff[64], m1[8][8], m2[8][8], m3[8][8], iSumHad = 0;

                    int x_pos = x + blkx;
                    if (x_pos >= width)
                        break;

                    Pel* piOrg = &curY[y_pos * stride + x_pos];
                    for (k = 0; k < 64; k += 8)
                    {
                        diff[k + 0] = piOrg[0];
                        diff[k + 1] = piOrg[1];
                        diff[k + 2] = piOrg[2];
                        diff[k + 3] = piOrg[3];
                        diff[k + 4] = piOrg[4];
                        diff[k + 5] = piOrg[5];
                        diff[k + 6] = piOrg[6];
                        diff[k + 7] = piOrg[7];

                        piOrg += stride;
                    }

                    //horizontal
                    for (j = 0; j < 8; j++)
                    {
                        jj = j << 3;
                        m2[j][0] = diff[jj    ] + diff[jj + 4];
                        m2[j][1] = diff[jj + 1] + diff[jj + 5];
                        m2[j][2] = diff[jj + 2] + diff[jj + 6];
                        m2[j][3] = diff[jj + 3] + diff[jj + 7];
                        m2[j][4] = diff[jj    ] - diff[jj + 4];
                        m2[j][5] = diff[jj + 1] - diff[jj + 5];
                        m2[j][6] = diff[jj + 2] - diff[jj + 6];
                        m2[j][7] = diff[jj + 3] - diff[jj + 7];

                        m1[j][0] = m2[j][0] + m2[j][2];
                        m1[j][1] = m2[j][1] + m2[j][3];
                        m1[j][2] = m2[j][0] - m2[j][2];
                        m1[j][3] = m2[j][1] - m2[j][3];
                        m1[j][4] = m2[j][4] + m2[j][6];
                        m1[j][5] = m2[j][5] + m2[j][7];
                        m1[j][6] = m2[j][4] - m2[j][6];
                        m1[j][7] = m2[j][5] - m2[j][7];

                        m2[j][0] = m1[j][0] + m1[j][1];
                        m2[j][1] = m1[j][0] - m1[j][1];
                        m2[j][2] = m1[j][2] + m1[j][3];
                        m2[j][3] = m1[j][2] - m1[j][3];
                        m2[j][4] = m1[j][4] + m1[j][5];
                        m2[j][5] = m1[j][4] - m1[j][5];
                        m2[j][6] = m1[j][6] + m1[j][7];
                        m2[j][7] = m1[j][6] - m1[j][7];
                    }

                    //vertical
                    for (i = 0; i < 8; i++)
                    {
                        m3[0][i] = m2[0][i] + m2[4][i];
                        m3[1][i] = m2[1][i] + m2[5][i];
                        m3[2][i] = m2[2][i] + m2[6][i];
                        m3[3][i] = m2[3][i] + m2[7][i];
                        m3[4][i] = m2[0][i] - m2[4][i];
                        m3[5][i] = m2[1][i] - m2[5][i];
                        m3[6][i] = m2[2][i] - m2[6][i];
                        m3[7][i] = m2[3][i] - m2[7][i];

                        m1[0][i] = m3[0][i] + m3[2][i];
                        m1[1][i] = m3[1][i] + m3[3][i];
                        m1[2][i] = m3[0][i] - m3[2][i];
                        m1[3][i] = m3[1][i] - m3[3][i];
                        m1[4][i] = m3[4][i] + m3[6][i];
                        m1[5][i] = m3[5][i] + m3[7][i];
                        m1[6][i] = m3[4][i] - m3[6][i];
                        m1[7][i] = m3[5][i] - m3[7][i];

                        m2[0][i] = m1[0][i] + m1[1][i];
                        m2[1][i] = m1[0][i] - m1[1][i];
                        m2[2][i] = m1[2][i] + m1[3][i];
                        m2[3][i] = m1[2][i] - m1[3][i];
                        m2[4][i] = m1[4][i] + m1[5][i];
                        m2[5][i] = m1[4][i] - m1[5][i];
                        m2[6][i] = m1[6][i] + m1[7][i];
                        m2[7][i] = m1[6][i] - m1[7][i];
                    }

                    for (i = 0; i < 8; i++)
                    {
                        for (j = 0; j < 8; j++)
                        {
                            iSumHad += abs(m2[i][j]);
                        }
                    }
                    iSumHad -= abs(m2[0][0]);
                    iSumHad = (iSumHad + 2) >> 2;
                    blk8_avg += iSumHad;
                    blk8_idx++;
                }
            }
            blk8_avg = (blk8_avg + blk8_idx - 1) / blk8_idx;
            encTopStat.CtuText[frame_sw][ctu_idx++].m_ctuHAD = blk8_avg;
            tmp += blk8_avg;
        }
    }
    encTopStat.FrameStat[frameno].m_avgHAD = (tmp + ctu_idx - 1) / ctu_idx;
}

#define ROBERTS_ED   0
void PreAnalyze_GRD(TComPic* curPic)
{
#define MATRIX_SZ       2
    TComPicYuv* pcPicYuv = curPic->getPicYuvOrg();
    const Pel* curY = pcPicYuv->getAddr(COMPONENT_Y);
    const int width = pcPicYuv->getWidth(COMPONENT_Y);
    const int height = pcPicYuv->getHeight(COMPONENT_Y);
    const int stride = pcPicYuv->getStride(COMPONENT_Y);
    const int ctb_size = curPic->getSlice(0)->getSPS()->getMaxCUWidth();
    const int frameno = encTopStat.m_encFrameNum;
    const int frame_sw = frameno & 0x1;

    int x, y, bx, by;
    int ctu_idx = 0;
    int tmp = 0;

    for (y = 0; y < height; y += ctb_size)
    {
        for (x = 0; x < width; x += ctb_size)
        {
            int x_size = ctb_size;
            int y_size = ctb_size;
            Int64 uiSum = 0;

            for (by = 0; by < ctb_size; by++)
            {
                int y_pos = y + by;
                if (y_pos >= (height - MATRIX_SZ + 1))
                {
                    y_size = by;
                    break;
                }
                for (bx = 0; bx < ctb_size; bx++)
                {
                    int x_pos = x + bx;
                    if (x_pos >= (width - MATRIX_SZ + 1))
                    {
                        x_size = bx;
                        break;
                    }
                    const Pel *pBlkY = &curY[y_pos * stride + x_pos];
#if ROBERTS_ED
                    int pix[4];

                    pix[0] = pBlkY[0];
                    pix[1] = pBlkY[1];
                    pix[2] = pBlkY[0 + stride];
                    pix[3] = pBlkY[1 + stride];

                    // [1  0]    [0  1]
                    // [0 -1]    [-1 0]
                    const Int64 grad0 = pix[0] - pix[3];
                    const Int64 grad1 = pix[1] - pix[2];

                    uiSum += grad0 * grad0 + grad1 * grad1;
#else
                    uiSum += abs(pBlkY[0] - pBlkY[1]) + abs(pBlkY[0] - pBlkY[0 + stride]);
#endif
                }
            }

#if ROBERTS_ED
            const int numPixInCTU = y_size * x_size;
            const Int64 dAverage = Int64(uiSum + numPixInCTU - 1) / numPixInCTU;
            encTopStat.CtuText[frame_sw][ctu_idx].m_ctuGRD = Int(dAverage);
            tmp += Int(dAverage);
            ctu_idx++;
#else
            encTopStat.CtuText[frame_sw][ctu_idx].m_ctuGRD = Int(uiSum);
            tmp += Int(uiSum);
            ctu_idx++;
#endif
        }
    }
    encTopStat.FrameStat[frameno].m_avgGRD = (tmp + ctu_idx - 1) / ctu_idx;
}

void PreAnalyze_SAD(TComPic* curPic, TComPic* refPic, int* diffPic)
{
    TComPicYuv* pcCurYuv = curPic->getPicYuvOrg();
    //TComPicYuv* pcRefYuv = refPic->getPicYuvRec();  // use reconstructed
    TComPicYuv* pcRefYuv = refPic->getPicYuvOrg();  // use source YUV
    const int width = pcCurYuv->getWidth(COMPONENT_Y);
    const int height = pcCurYuv->getHeight(COMPONENT_Y);
    const int stride = pcCurYuv->getStride(COMPONENT_Y);
    const int ctb_size = curPic->getSlice(0)->getSPS()->getMaxCUWidth();
    const Pel* curY = pcCurYuv->getAddr(COMPONENT_Y);
    const Pel* refY = pcRefYuv->getAddr(COMPONENT_Y);
    const int frameno = encTopStat.m_encFrameNum;
    const int frame_sw = frameno & 0x1;

    assert(pcCurYuv->getWidth(COMPONENT_Y) == pcRefYuv->getWidth(COMPONENT_Y));
    assert(pcCurYuv->getHeight(COMPONENT_Y) == pcRefYuv->getHeight(COMPONENT_Y));
    assert(pcCurYuv->getStride(COMPONENT_Y) == pcRefYuv->getStride(COMPONENT_Y));

    int x, y, bx, by;
    int ctu_idx = 0;
    Int64 tmp = 0;

    for (y = 0; y < height; y += ctb_size)
    {
        for (x = 0; x < width; x += ctb_size)
        {
            int x_size = ctb_size;
            int y_size = ctb_size;
            int uiSum = 0;
            int uiSumC = 0;

#if 1   // Find Spatial DC
            int DC_val = 0;
            int pix_cnt = 0;
            if (x > 1)
            {
                for (by = 0; by < ctb_size; by++)
                {
                    const int y_pos = y + by;
                    if (y_pos >= height)
                        break;
                    const int pos = y_pos * stride + (x - 1);
                    DC_val += curY[pos];
                    pix_cnt++;
                }
            }
            else
            {
                DC_val += 128*ctb_size;
                pix_cnt += ctb_size;
            }

            if (y > 1)
            {
                for (bx = 0; bx < ctb_size; bx++)
                {
                    const int x_pos = x + bx;
                    if (x_pos >= width)
                        break;
                    const int pos = (y - 1) * stride + x_pos;
                    DC_val += curY[pos];
                    pix_cnt++;
                }
            }
            else
            {
                DC_val += 128 * ctb_size;
                pix_cnt += ctb_size;
            }
            DC_val = (DC_val + (pix_cnt - 1)) / pix_cnt;
#endif

            for (by = 0; by < ctb_size; by++)
            {
                int y_pos = y + by;
                if (y_pos >= height)
                {
                    y_size = by;
                    break;
                }
                for (bx = 0; bx < ctb_size; bx++)
                {
                    int x_pos = x + bx;
                    if (x_pos >= width)
                    {
                        x_size = bx;
                        break;
                    }
                    const int pos = y_pos * stride + x_pos;
                    const int posD = y_pos * width + x_pos;
                    const Pel pCurY = curY[pos];
                    const Pel pRefY = refY[pos];
                    diffPic[posD] = pCurY - pRefY;
                    uiSum += abs(diffPic[posD]);
                    //uiSum2 += diffPic[posD] * diffPic[posD];

                    uiSumC += abs(pCurY - DC_val);

                    //uiSumS += pCurY;
                    //uiSumS2 += pCurY * pCurY;
                }
            }

#if 0   // Est Skip Block and Bits
            uiSum = max(uiSum, 1);  // prevent to divide 0
            // Estimate Skip Blocks
            int a0;
            const int refCtuSAD = encTopStat.CtuText[ref_sw][ctu_idx].m_ctuSAD;
            const int refSkips = encTopStat.CtuText[ref_sw][ctu_idx].m_CtuType[0];
            const int refPs = encTopStat.CtuText[ref_sw][ctu_idx].m_CtuType[1];
            if (refCtuSAD > uiSum) // curSAD is lower, guess skips blks will be more
                if(refSkips == 0)
                    a0 = int(refPs * refCtuSAD / uiSum / 5);
                else
                    a0 = refSkips + int(refPs * refCtuSAD / uiSum);	// there might be intra?
            else // curSAD is higher, guess skips blks will be less
            {
                if ((uiSum < 2048) || ((uiSum - refCtuSAD) < 1024))
                    a0 = (refSkips == 0) ? refPs : refSkips;
                else if (uiSum < 4096)
                    a0 = int(refSkips * (2 - uiSum/2048));
                else
                    a0 = int(refSkips * refCtuSAD / uiSum);
            }
            encTopStat.CtuText[frame_sw][ctu_idx].m_estSkip = min(a0, total_8x8);
#endif
            uiSum = min(uiSum, uiSumC);
            encTopStat.CtuText[frame_sw][ctu_idx].m_ctuSAD = uiSum;

            /*const int numPixInCTU = y_size * x_size;
            const Int64 dAverage = Int64(uiSum + numPixInCTU - 1) / numPixInCTU;
            const Int64 dAverageS = Int64(uiSumS + numPixInCTU - 1) / numPixInCTU;
            //const Int64 dVariance = Int64(uiSum2 + numPixInCTU - 1) / numPixInCTU - dAverage * dAverage + 1;
            const Int64 dVarianceS = abs(uiSumS2 - uiSumS * dAverageS) + 1;
            const Int64 dVariance = min(abs(uiSum2 - uiSum * dAverage) + 1, dVarianceS);

            encTopStat.CtuText[frame_sw][ctu_idx].m_ctuDVAR = Int(dVariance);*/
                        
            tmp += uiSum;
            ctu_idx++;
        }
    }
    assert(ctu_idx == refPic->getPicSym()->getNumberOfCtusInFrame());
    encTopStat.FrameStat[frameno].m_TotalSAD = tmp;
    encTopStat.FrameStat[frameno].m_avgSAD = Int((tmp + ctu_idx - 1) / ctu_idx);
    //encTopStat.FrameStat[frameno].m_avgDVAR = Int((tmp3 + ctu_idx - 1) / ctu_idx);
}

void PreAnalyze_DHAD(TComPic* curPic, int* diffPic)
{
    TComPicYuv* pcPicYuv = curPic->getPicYuvOrg();
    //Pel* curY = pcPicYuv->getAddr(COMPONENT_Y);
    const int width = pcPicYuv->getWidth(COMPONENT_Y);
    const int height = pcPicYuv->getHeight(COMPONENT_Y);
    //const int stride = pcPicYuv->getStride(COMPONENT_Y);
    const int ctb_size = curPic->getSlice(0)->getSPS()->getMaxCUWidth();
    const int blk_num = ctb_size >> 3;
    const int frameno = encTopStat.m_encFrameNum;
    const int frame_sw = frameno & 0x1;

    int x, y, blkx, blky;
    int ctu_idx = 0;
    int tmp = 0;

    for (y = 0; y < height; y += ctb_size)
    {
        for (x = 0; x < width; x += ctb_size)
        {
            int blk8_idx = 0;
            int blk8_avg = 0;
            for (blky = 0; blky < ctb_size; blky += 8)  // frame size must align to 8
            {
                int y_pos = y + blky;
                if (y_pos >= height)
                    break;

                for (blkx = 0; blkx < ctb_size; blkx += 8)
                {
                    int k, i, j, jj;
                    int diff[64], m1[8][8], m2[8][8], m3[8][8], iSumHad = 0;

                    int x_pos = x + blkx;
                    if (x_pos >= width)
                        break;

                    int* piOrg = &diffPic[y_pos * width + x_pos];
                    for (k = 0; k < 64; k += 8)
                    {
                        diff[k + 0] = abs(piOrg[0]);
                        diff[k + 1] = abs(piOrg[1]);
                        diff[k + 2] = abs(piOrg[2]);
                        diff[k + 3] = abs(piOrg[3]);
                        diff[k + 4] = abs(piOrg[4]);
                        diff[k + 5] = abs(piOrg[5]);
                        diff[k + 6] = abs(piOrg[6]);
                        diff[k + 7] = abs(piOrg[7]);

                        piOrg += width;
                    }

                    //horizontal
                    for (j = 0; j < 8; j++)
                    {
                        jj = j << 3;
                        m2[j][0] = diff[jj] + diff[jj + 4];
                        m2[j][1] = diff[jj + 1] + diff[jj + 5];
                        m2[j][2] = diff[jj + 2] + diff[jj + 6];
                        m2[j][3] = diff[jj + 3] + diff[jj + 7];
                        m2[j][4] = diff[jj] - diff[jj + 4];
                        m2[j][5] = diff[jj + 1] - diff[jj + 5];
                        m2[j][6] = diff[jj + 2] - diff[jj + 6];
                        m2[j][7] = diff[jj + 3] - diff[jj + 7];

                        m1[j][0] = m2[j][0] + m2[j][2];
                        m1[j][1] = m2[j][1] + m2[j][3];
                        m1[j][2] = m2[j][0] - m2[j][2];
                        m1[j][3] = m2[j][1] - m2[j][3];
                        m1[j][4] = m2[j][4] + m2[j][6];
                        m1[j][5] = m2[j][5] + m2[j][7];
                        m1[j][6] = m2[j][4] - m2[j][6];
                        m1[j][7] = m2[j][5] - m2[j][7];

                        m2[j][0] = m1[j][0] + m1[j][1];
                        m2[j][1] = m1[j][0] - m1[j][1];
                        m2[j][2] = m1[j][2] + m1[j][3];
                        m2[j][3] = m1[j][2] - m1[j][3];
                        m2[j][4] = m1[j][4] + m1[j][5];
                        m2[j][5] = m1[j][4] - m1[j][5];
                        m2[j][6] = m1[j][6] + m1[j][7];
                        m2[j][7] = m1[j][6] - m1[j][7];
                    }

                    //vertical
                    for (i = 0; i < 8; i++)
                    {
                        m3[0][i] = m2[0][i] + m2[4][i];
                        m3[1][i] = m2[1][i] + m2[5][i];
                        m3[2][i] = m2[2][i] + m2[6][i];
                        m3[3][i] = m2[3][i] + m2[7][i];
                        m3[4][i] = m2[0][i] - m2[4][i];
                        m3[5][i] = m2[1][i] - m2[5][i];
                        m3[6][i] = m2[2][i] - m2[6][i];
                        m3[7][i] = m2[3][i] - m2[7][i];

                        m1[0][i] = m3[0][i] + m3[2][i];
                        m1[1][i] = m3[1][i] + m3[3][i];
                        m1[2][i] = m3[0][i] - m3[2][i];
                        m1[3][i] = m3[1][i] - m3[3][i];
                        m1[4][i] = m3[4][i] + m3[6][i];
                        m1[5][i] = m3[5][i] + m3[7][i];
                        m1[6][i] = m3[4][i] - m3[6][i];
                        m1[7][i] = m3[5][i] - m3[7][i];

                        m2[0][i] = m1[0][i] + m1[1][i];
                        m2[1][i] = m1[0][i] - m1[1][i];
                        m2[2][i] = m1[2][i] + m1[3][i];
                        m2[3][i] = m1[2][i] - m1[3][i];
                        m2[4][i] = m1[4][i] + m1[5][i];
                        m2[5][i] = m1[4][i] - m1[5][i];
                        m2[6][i] = m1[6][i] + m1[7][i];
                        m2[7][i] = m1[6][i] - m1[7][i];
                    }

                    for (i = 0; i < 8; i++)
                    {
                        for (j = 0; j < 8; j++)
                        {
                            iSumHad += abs(m2[i][j]);
                        }
                    }
                    iSumHad -= abs(m2[0][0]);
                    iSumHad = (iSumHad + 2) >> 2;
                    blk8_avg += iSumHad;
                    blk8_idx++;
                }
            }
            blk8_avg = (blk8_avg + blk8_idx - 1) / blk8_idx;
            encTopStat.CtuText[frame_sw][ctu_idx].m_ctuDHAD = blk8_avg;
            tmp += blk8_avg;
            ctu_idx++;
        }
    }
    encTopStat.FrameStat[frameno].m_avgDHAD = (tmp + ctu_idx - 1) / ctu_idx;
}

void PreAnalyze_DGR(TComPic* curPic, int* diffPic)
{
#define MATRIX_SZ       2
    TComPicYuv* pcPicYuv = curPic->getPicYuvOrg();
    //TComPicYuv* pcRefYuv = refPic->getPicYuvOrg();
    const int width = pcPicYuv->getWidth(COMPONENT_Y);
    const int height = pcPicYuv->getHeight(COMPONENT_Y);
    const int stride = pcPicYuv->getStride(COMPONENT_Y);
    const int ctb_size = curPic->getSlice(0)->getSPS()->getMaxCUWidth();
    const Pel* curY = pcPicYuv->getAddr(COMPONENT_Y);
    //const Pel* refY = pcRefYuv->getAddr(COMPONENT_Y);
    const int frameno = encTopStat.m_encFrameNum;
    const int frame_sw = frameno & 0x1;

    int x, y, bx, by;
    int ctu_idx = 0;
    Int64 tmp = 0;

    for (y = 0; y < height; y += ctb_size)
    {
        for (x = 0; x < width; x += ctb_size)
        {
            int x_size = ctb_size;
            int y_size = ctb_size;
            int uiSum = 0;
            int uiSum1 = 0;

            for (by = 0; by < ctb_size; by++)
            {
                int y_pos = y + by;
                if (y_pos >= (height - MATRIX_SZ + 1))
                {
                    y_size = by;
                    break;
                }
                for (bx = 0; bx < ctb_size; bx++)
                {
                    int x_pos = x + bx;
                    if (x_pos >= (width - MATRIX_SZ + 1))
                    {
                        x_size = bx;
                        break;
                    }
                    const int* pBlkY = &diffPic[y_pos * width + x_pos];
#if ROBERTS_ED
                    int pix[4];

                    pix[0] = pBlkY[0];
                    pix[1] = pBlkY[1];
                    pix[2] = pBlkY[0 + width];
                    pix[3] = pBlkY[1 + width];

                    // [1  0]    [0  1]
                    // [0 -1]    [-1 0]
                    const Int64 grad0 = pix[0] - pix[3];
                    const Int64 grad1 = pix[1] - pix[2];

                    uiSum += grad0 * grad0 + grad1 * grad1;
#else
                    const Pel* pSrcY = &curY[y_pos * stride + x_pos];
                    //const Pel* pRefY = &refY[y_pos * stride + x_pos];
                    uiSum += abs(pBlkY[0] - pBlkY[1]) + abs(pBlkY[0] - pBlkY[0 + width]);   // temporal
                    //uiSum  += abs(pRefY[0] - pRefY[1]) + abs(pRefY[0] - pRefY[0 + stride]);   // ref, spatial
                    uiSum1 += abs(pSrcY[0] - pSrcY[1]) + abs(pSrcY[0] - pSrcY[0 + stride]); // cur, spatial
#endif
                }
            }

#if ROBERTS_ED
            //const int numPixInCTU = y_size * x_size;
            //const int dAverage = (uiSum + numPixInCTU - 1) / numPixInCTU;
            //encTopStat.CtuText[frame_sw][ctu_idx].m_ctuDGR = dAverage;
            //tmp += dAverage;
            encTopStat.CtuText[frame_sw][ctu_idx].m_ctuDGR = uiSum;
            encTopStat.FrameStat[frameno].m_MotionDGR[encTopStat.CtuText[frame_sw][ctu_idx].m_motion] += uiSum;
            ctu_idx++;
            tmp += uiSum;
#else
            int minGRD = min(uiSum, uiSum1);
            //int minGRD = min(abs(uiSum-uiSum1), uiSum1);
            encTopStat.CtuText[frame_sw][ctu_idx].m_ctuDGR = minGRD;
            tmp += minGRD;
            ctu_idx++;
#endif
        }
    }
    encTopStat.FrameStat[frameno].m_avgDGR = int(tmp + ctu_idx - 1) / ctu_idx;
}

extern double HAD_ALPHA_LIN[52];
extern double HAD_BETA_LIN[52];
void Enc_Init_Text(int width, int height, int ctu_size)
{
    const int frame_width_ctu = (width + ctu_size - 1) / ctu_size;
    const int frame_height_ctu = (height + ctu_size - 1) / ctu_size;
    const int ctu_num = frame_width_ctu * frame_height_ctu;

    encTopStat.m_encFrameNum = 0;
    encTopStat.CtuText[0] = new EncCtuText[ctu_num];
    encTopStat.CtuText[1] = new EncCtuText[ctu_num];
    memset(encTopStat.CtuText[0], 0, sizeof(EncCtuText) * ctu_num);
    memset(encTopStat.CtuText[1], 0, sizeof(EncCtuText) * ctu_num);

    // Initialize RC parameters
    memcpy(encTopStat.HAD_ALPHA, HAD_ALPHA_LIN, sizeof(HAD_ALPHA_LIN));
    memcpy(encTopStat.HAD_BETA,  HAD_BETA_LIN, sizeof(HAD_BETA_LIN));

    system("mkdir ctu_info");
}

void Enc_DeInit_Text()
{
    delete[] encTopStat.CtuText[0];
    delete[] encTopStat.CtuText[1];
}

#define NORM_ONE    1   // 1: normalize to 1 and clip, plus a shift, 0: just scale to clip range
void PreAnalyze_CTU_Weight(TComPic* curPic)
{
    TComPicYuv* pcPicYuv = curPic->getPicYuvOrg();
    const int width = pcPicYuv->getWidth(COMPONENT_Y);
    const int height = pcPicYuv->getHeight(COMPONENT_Y);
    const int stride = pcPicYuv->getStride(COMPONENT_Y);
    const int ctb_size = curPic->getSlice(0)->getSPS()->getMaxCUWidth();
    const Pel* curY = pcPicYuv->getAddr(COMPONENT_Y);
    const int frameno = encTopStat.m_encFrameNum;
    const int frame_sw = frameno & 0x1;

    int x, y, bx, by, y_pos, x_pos, round;
    int ctu_idx = 0;

    double uiSumW = 0;
    double uiMaxV = 0;
    double uiMinV = 0xFFFFFF;
    UInt64 uiSumY = 0;
    UInt64 uiSumUV = 0;
    double apic = 16 * pow(3840 * 2160 / width / height, 0.25);
    Pel p0, p1, p2, p3, p4, p5, p6, p7, p8, plow, phigh;
    // Low-pass Luma Pic
    for (y = 0; y < height; y += ctb_size)
    {
        for (x = 0; x < width; x += ctb_size)
        {
            int x_size = ctb_size;
            int y_size = ctb_size;
            int uiSumH = 0;
            int uiSumLuma = 0;
            int uiSumLuma2 = 0;

            for (by = 0; by < ctb_size; by++)
            {
                if ((y + by) >= height)
                {
                    y_size = by;
                    break;
                }
                for (bx = 0; bx < ctb_size; bx++)
                {
                    if ((x + bx) >= width)
                    {
                        x_size = bx;
                        break;
                    }
#if CTU_WEIGHT==1
                    y_pos = Clip3(0, height - 1, y + by - 1);
                    x_pos = Clip3(0, width - 1, x + bx - 1);
                    p0 = curY[y_pos * stride + x_pos];
                    x_pos = Clip3(0, width - 1, x + bx + 0);
                    p1 = curY[y_pos * stride + x_pos];
                    x_pos = Clip3(0, width - 1, x + bx + 1);
                    p2 = curY[y_pos * stride + x_pos];

                    y_pos = Clip3(0, height - 1, y + by + 0);
                    x_pos = Clip3(0, width - 1, x + bx - 1);
                    p3 = curY[y_pos * stride + x_pos];
                    x_pos = Clip3(0, width - 1, x + bx + 0);
                    p4 = curY[y_pos * stride + x_pos];
                    x_pos = Clip3(0, width - 1, x + bx + 1);
                    p5 = curY[y_pos * stride + x_pos];

                    y_pos = Clip3(0, height - 1, y + by + 1);
                    x_pos = Clip3(0, width - 1, x + bx - 1);
                    p6 = curY[y_pos * stride + x_pos];
                    x_pos = Clip3(0, width - 1, x + bx + 0);
                    p7 = curY[y_pos * stride + x_pos];
                    x_pos = Clip3(0, width - 1, x + bx + 1);
                    p8 = curY[y_pos * stride + x_pos];

                    plow = (p0 + 2 * p1 + p2 + 2 * p3 + 4 * p4 + 2 * p5 + p6 + 2 * p7 + p8) >> 2;   // use 4 ?
                    phigh = abs(4 * p4 - plow);
                    uiSumH += phigh;
                    uiSumLuma += p4;    // for Luma weight
#else
                    y_pos = Clip3(0, height - 1, y + by + 0);
                    x_pos = Clip3(0, width - 1, x + bx + 0);
                    p4 = curY[y_pos * stride + x_pos];

                    uiSumLuma += p4;
                    uiSumLuma2 += p4 * p4;
#endif
                }
            }

            const int numPixInCTU = y_size * x_size;
#if CTU_WEIGHT==1
            // CTU_WEIGHT = 1
            uiSumY += uiSumH;
            double hk = (double)uiSumH / numPixInCTU;
            double ak = max(4.0, hk);    // 4 is lower sensitivity limit of HVS
            double wk = apic / ak;

            /*int cur_luma = (uiSumLuma + numPixInCTU - 1) / numPixInCTU;
            encTopStat.CtuText[frame_sw][ctu_idx].m_ctuAvgLuma = cur_luma;
            if (cur_luma < 75)
                encTopStat.CtuText[frame_sw][ctu_idx].luma_weight = (double)75 / (150 - cur_luma);
            else if (cur_luma <= 125)
                encTopStat.CtuText[frame_sw][ctu_idx].luma_weight = 1;
            else
                encTopStat.CtuText[frame_sw][ctu_idx].luma_weight = (double)92 / (cur_luma - 32);
            wk *= encTopStat.CtuText[frame_sw][ctu_idx].luma_weight;*/

            encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight = wk;
            uiSumW += wk;
            if (wk > uiMaxV)
                uiMaxV = wk;
            else if (wk < uiMinV)
                uiMinV = wk;
#else
            // CTU_WEIGHT = 2, 3
            const int dAverage = (uiSumLuma + numPixInCTU - 1) / numPixInCTU;
            const int dVariance = abs(uiSumLuma2 - uiSumLuma * dAverage) + 1;
            const int wki = (dVariance + numPixInCTU - 1) / numPixInCTU;
            encTopStat.CtuText[frame_sw][ctu_idx].m_ctuAvgLuma = dAverage;
            encTopStat.CtuText[frame_sw][ctu_idx].m_ctuAvgVAR = wki;
            uiSumW += wki;
            if (wki > uiMaxV)
                uiMaxV = wki;
            else if (wki < uiMinV)
                uiMinV = wki;
#endif
            ctu_idx++;
        }
    }

#if CTU_WEIGHT==1
    // Normalize to 1
    // clipping the range equals to limit QP
    // QP [-6,-5,-4,-3,-2,-1, 0] = [3.93, 3.10, 2.44, 1.92, 1.51, 1.19, 1.0]
    // QP [ 6, 5, 4, 3, 2, 1, 0] = [0.25, 0.32, 0.41, 0.52, 0.66, 0.84, 1.0]
    double QP_NEG[7] = { 4.17, 3.28, 2.59, 2.04, 1.60, 1.26, 1.00 };
    double QP_POS[7] = { 0.24, 0.30, 0.38, 0.49, 0.62, 0.79, 1.00 };
#if 1   // divide
    double QPMax = QP_NEG[6 - 4];   // QP - 4
    double QPMin = QP_POS[6 - 2];   // QP + 2
#else   // miltiply
    double QPMax = QP_NEG[6 - 2];   // QP + 2
    double QPMin = QP_POS[6 - 4];   // QP - 4
#endif

    const int ctu_num = curPic->getNumberOfCtusInFrame();
#if NORM_ONE
    double adjust = ctu_num / uiSumW;  // rescale to close to 1
    double shift = 0;
    //double adjust = encTopCRC.CRCCtuLambdaScale * ctu_idx / uiSumW;
    double uiSumW2 = 0;
    for (ctu_idx = 0; ctu_idx < ctu_num; ctu_idx++)
    {
        double wk = Clip3(QPMin, QPMax, encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight * adjust);    // Pic QP +-X
        //double wk = encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight * adjust;    // Pic QP +-X
        //encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight = wk;
        uiSumW2 += wk;

        #if 1 // divide , 1/log(1.0613235) = 16.802
        encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight = -1 * int(16.802 * log(wk) + 0.5);
        #else // multiply
        encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight = int(16.802 * log(wk) + 0.5);
        #endif
    }
#else
    double adjust = (QPMax - QPMin) / (uiMaxV - uiMinV); // rescale to close to QP-delta range
    double shift = QPMax - uiMaxV * adjust;
    double uiSumW2 = 0;
    for (ctu_idx = 0; ctu_idx < ctu_num; ctu_idx++)
    {
        double wk = encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight * adjust + shift;    // Pic QP +-X
        //encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight = wk;
        uiSumW2 += wk;

        #if 1 // divide , 1/log(1.0613235) = 16.802
        encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight = -1 * int(16.802 * log(wk) + 0.5);
        #else // multiply
        encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight = int(16.802 * log(wk) + 0.5);
        #endif
    }
#endif
    // check average weight close to 1
    double is_1 = uiSumW2 / ctu_num;
    printf("RC: Highpass average weight = %f, adjust = %f, shift = %f\n", is_1, adjust, shift);

    // Low-pass Chroma Pic
    const int width_c = pcPicYuv->getWidth(COMPONENT_Cb);
    const int height_c = pcPicYuv->getHeight(COMPONENT_Cb);
    const int stride_c = pcPicYuv->getStride(COMPONENT_Cb);
    const int ctb_size_c = ctb_size / 2;
    for (round = 0; round < 2; round++)
    {
        const Pel* curUV = (round == 0) ? pcPicYuv->getAddr(COMPONENT_Cb) : pcPicYuv->getAddr(COMPONENT_Cr);
        for (y = 0; y < height_c; y++)
        {
            for (x = 0; x < width_c; x++)
            {
                y_pos = Clip3(0, height_c - 1, y - 1);
                x_pos = Clip3(0, width_c - 1, x - 1);
                p0 = curUV[y_pos * stride_c + x_pos];
                x_pos = Clip3(0, width_c - 1, x + 0);
                p1 = curUV[y_pos * stride_c + x_pos];
                x_pos = Clip3(0, width_c - 1, x + 1);
                p2 = curUV[y_pos * stride_c + x_pos];

                y_pos = Clip3(0, height_c - 1, y + 0);
                x_pos = Clip3(0, width_c - 1, x - 1);
                p3 = curUV[y_pos * stride_c + x_pos];
                x_pos = Clip3(0, width_c - 1, x + 0);
                p4 = curUV[y_pos * stride_c + x_pos];
                x_pos = Clip3(0, width_c - 1, x + 1);
                p5 = curUV[y_pos * stride_c + x_pos];

                y_pos = Clip3(0, height_c - 1, y + 1);
                x_pos = Clip3(0, width_c - 1, x - 1);
                p6 = curUV[y_pos * stride_c + x_pos];
                x_pos = Clip3(0, width_c - 1, x + 0);
                p7 = curUV[y_pos * stride_c + x_pos];
                x_pos = Clip3(0, width_c - 1, x + 1);
                p8 = curUV[y_pos * stride_c + x_pos];

                plow = (p0 + 2 * p1 + p2 + 2 * p3 + 4 * p4 + 2 * p5 + p6 + 2 * p7 + p8) >> 2;   // use 4 ?
                phigh = abs(4*p4 - plow);
                uiSumUV += phigh;
            }
        }

        if (8 * uiSumUV > uiSumY)
        {
            encTopCRC.CRCQPC_delta[round] = min(4, int(3 * log2(8 * uiSumUV / uiSumY) + 0.5));  // 3 * log2(wk)
            encTopStat.FrameStat[frameno].blk_weight_c[round] = 8 * (double)uiSumUV / uiSumY;
            printf("RC: Chroma Weight[%d] = %f, delta = %d\n", round, encTopStat.FrameStat[frameno].blk_weight_c[round], encTopCRC.CRCQPC_delta[round]);
        }
        else
        {
            encTopCRC.CRCQPC_delta[round] = 0;
            encTopStat.FrameStat[frameno].blk_weight_c[round] = 0;
        }
    }
#elif CTU_WEIGHT == 2
    // Calculate CTU lambda weight
    double weightSum = 0;
    const int ctu_num = curPic->getNumberOfCtusInFrame();
    const UInt64 total_sad = encTopStat.FrameStat[frameno].m_TotalSAD;  // make sure this available
    int scaleMotion = 1024;

    uiMaxV = 0;
    uiMinV = 0xFFFFFF;
    for (ctu_idx = 0; ctu_idx < ctu_num; ctu_idx++)
    {
        int cur_luma = encTopStat.CtuText[frame_sw][ctu_idx].m_ctuAvgLuma;
        if (cur_luma < 75)
            encTopStat.CtuText[frame_sw][ctu_idx].luma_weight = (double)75 / (150 - cur_luma);
        else if (cur_luma <= 125)
            encTopStat.CtuText[frame_sw][ctu_idx].luma_weight = 1;
        else
            encTopStat.CtuText[frame_sw][ctu_idx].luma_weight = (double)92 / (cur_luma - 32);

        int cur_sad = encTopStat.CtuText[frame_sw][ctu_idx].m_ctuSAD;
        //encTopStat.CtuText[frame_sw][ctu_idx].motion_weight = (double)cur_sad * scaleMotion / (cur_sad + total_sad);
        //encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight = encTopStat.CtuText[frame_sw][ctu_idx].luma_weight * encTopStat.CtuText[frame_sw][ctu_idx].motion_weight;
        encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight = encTopStat.CtuText[frame_sw][ctu_idx].luma_weight * cur_sad;
        weightSum += encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight;

        if (encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight > uiMaxV)
            uiMaxV = encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight;
        else if (encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight < uiMinV)
            uiMinV = encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight;
    }

    /*double weightSum2 = 0;
    for (ctu_idx = 0; ctu_idx < ctu_num; ctu_idx++)
    {
        Double CtuWeight = encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight / weightSum;
        //Double K = 1.2 * CtuWeight + 0.3;
        Double K = CtuWeight;
        encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight = K;
        weightSum2 += K;
    }*/

    // Normalize to 1
    // clipping the range equals to limit QP
    // QP [ 6, 5, 4, 3, 2, 1, 0] = [3.93, 3.10, 2.44, 1.92, 1.51, 1.19, 1.0]
    // QP [-6,-5,-4,-3,-2,-1, 0] = [0.25, 0.32, 0.41, 0.52, 0.66, 0.84, 1.0]
    double QP_POS[7] = { 4.17, 3.28, 2.59, 2.04, 1.60, 1.26, 1.00 };
    double QP_NEG[7] = { 0.24, 0.30, 0.38, 0.49, 0.62, 0.79, 1.00 };
    double QPMax = QP_POS[6 - 2];   // QP + 2
    double QPMin = QP_NEG[6 - 2];   // QP - 2
    double uiSumW2 = 0;
    //double adjust1 = ctu_num / weightSum2 / encTopCRC.CRCCtuLambdaScale;
    double adjust2 = (QPMax - QPMin) / (uiMaxV - uiMinV); // rescale to close to QP-delta range
    double adjust1 = ctu_num / weightSum;  // rescale to close to 1
    //double adjust3 = adjust2 / adjust1;
    double shift3 = QPMax - uiMaxV * adjust2;
    for (ctu_idx = 0; ctu_idx < ctu_num; ctu_idx++)
    {
        // double wk = Clip3(QPMin, QPMax, encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight * adjust1);
        //double wk = encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight * adjust1;
        double wk = encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight * adjust2 + shift3;    // Pic QP +-X
        //encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight = wk;
        uiSumW2 += wk;

        encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight = int(16.802 * log(wk) + 0.5);
    }

    // check average weight close to 1
    double is_1 = uiSumW2 / ctu_num;
    //printf("RC: Luma/Motion average weight = %f, adjust= %f\n", is_1, adjust1);
    printf("RC: Luma/Motion average weight = %f, adjust = %f, shift = %f\n", is_1, adjust2, shift3);
#elif CTU_WEIGHT == 3
    // Normalize to 1
    // clipping the range equals to limit QP
    // QP [ 6, 5, 4, 3, 2, 1, 0] = [3.93, 3.10, 2.44, 1.92, 1.51, 1.19, 1.0]
    // QP [-6,-5,-4,-3,-2,-1, 0] = [0.25, 0.32, 0.41, 0.52, 0.66, 0.84, 1.0]
    double QP_POS[7] = { 4.17, 3.28, 2.59, 2.04, 1.60, 1.26, 1.00 };
    double QP_NEG[7] = { 0.24, 0.30, 0.38, 0.49, 0.62, 0.79, 1.00 };
    double QPMax = QP_POS[6 - 2];   // QP + 2
    double QPMin = QP_NEG[6 - 2];   // QP - 2
    double uiSumW2 = 0;

    const int ctu_num = curPic->getNumberOfCtusInFrame();
    double adjust2 = (QPMax - QPMin) / (uiMaxV - uiMinV); // rescale to close to QP-delta range
    double adjust1 = ctu_num / uiSumW;  // rescale to close to 1
    //double adjust3 = adjust2 / adjust1;
    double shift3 = QPMax - uiMaxV * adjust2;
    for (ctu_idx = 0; ctu_idx < ctu_num; ctu_idx++)
    {
        //double wk = encTopStat.CtuText[frame_sw][ctu_idx].var_weight * adjust2 + shift3;    // Pic QP +-X
        double wk = encTopStat.CtuText[frame_sw][ctu_idx].m_ctuAvgVAR * adjust2 + shift3;    // Pic QP +-X
        //encTopStat.CtuText[frame_sw][ctu_idx].var_weight = wk;
        uiSumW2 += wk;

        encTopStat.CtuText[frame_sw][ctu_idx].var_weight = int(16.802 * log(wk) + 0.5);
    }

    // check average weight close to 1
    double is_1 = uiSumW2 / ctu_num;
    //printf("RC: Variance average weight = %f, adjust= %f\n", is_1, adjust1);
    printf("RC: Variance average weight = %f, adjust = %f, shift = %f\n", is_1, adjust2, shift3);
#endif
}


void Analyze_Pictures(TComPic * curPic, int intra_period)
{
    EncFrameStat framestat;
    clock_t iBeforeTime;
    Double m_EndTime;
    memset(&framestat, 0, sizeof(EncFrameStat));
    encTopStat.FrameStat.push_back(framestat);

    if (curPic->getSlice(0)->getSliceType() == I_SLICE)
    {
        iBeforeTime = clock();
        PreAnalyze_HAD(curPic);
        m_EndTime = (Double)(clock() - iBeforeTime) / CLOCKS_PER_SEC;
        //printf("Frame %d HAD Time %6.3f\n", encTopStat.m_encFrameNum, m_EndTime);

#if CTU_WEIGHT==2
        const int frameno = encTopStat.m_encFrameNum;
        const int frame_sw = frameno & 0x1;
        const int ctu_num = curPic->getNumberOfCtusInFrame();
        for (int ctu_idx = 0; ctu_idx < ctu_num; ctu_idx++)
        {
            encTopStat.CtuText[frame_sw][ctu_idx].lumo_weight = 1;
        }
#else
        /*const int frameno = encTopStat.m_encFrameNum;
        if (frameno == 0)
        {
            const int frame_sw = frameno & 0x1;
            const int ctu_num = curPic->getNumberOfCtusInFrame();
            for (int ctu_idx = 0; ctu_idx < ctu_num; ctu_idx++)
            {
                encTopStat.CtuText[frame_sw][ctu_idx].hipa_weight = 1;
            }
        }
        else */if (encTopCRC.CRCEn && encTopCRC.CRCCtuEn)
        {
            PreAnalyze_CTU_Weight(curPic);
        }
#endif
    }
    else //if (curPic->getSlice(0)->getSliceType() == P_SLICE)
    {
        const int width = curPic->getPicYuvOrg()->getWidth(COMPONENT_Y);
        const int height = curPic->getPicYuvOrg()->getHeight(COMPONENT_Y);
        int* diffPic = new int[width * height];

        TComPic* refPic = encTopStat.refPic;

        iBeforeTime = clock();
        PreAnalyze_SAD(curPic, refPic, diffPic);
        m_EndTime = (Double)(clock() - iBeforeTime) / CLOCKS_PER_SEC;
        //printf("Frame %d SAD Time %6.3f\n", encTopStat.m_encFrameNum, m_EndTime);

        /*iBeforeTime = clock();
        PreAnalyze_DHAD(curPic, diffPic);
        m_EndTime = (Double)(clock() - iBeforeTime) / CLOCKS_PER_SEC;
        printf("Frame %d DHAD Time %6.3f\n", encTopStat.m_encFrameNum, m_EndTime);*/
        
        iBeforeTime = clock();
        PreAnalyze_DGR(curPic, diffPic);
        m_EndTime = (Double)(clock() - iBeforeTime) / CLOCKS_PER_SEC;
        //printf("Frame %d DGR Time %6.3f\n", encTopStat.m_encFrameNum, m_EndTime);

        if ((encTopStat.m_encFrameNum % intra_period) == (intra_period-1))
        {
            iBeforeTime = clock();
            PreAnalyze_HAD(curPic);
            m_EndTime = (Double)(clock() - iBeforeTime) / CLOCKS_PER_SEC;
            //printf("Frame %d HAD Time %6.3f\n", encTopStat.m_encFrameNum, m_EndTime);
        }

        delete[] diffPic;
#if REDUCED_ENCODER_MEMORY && CUSTOM_RC 
        //if (!isField) // don't release the source data for field-coding because the fields are dealt with in pairs. // TODO: release source data for interlace simulations.
        {
            refPic->releaseEncoderSourceImageData();
        }
#endif
        if (encTopCRC.CRCEn && encTopCRC.CRCCtuEn)
        {
            PreAnalyze_CTU_Weight(curPic);
        }
    }

    encTopStat.refPic = curPic;    // record last YUV for I-P-P-P (only?)
}


void Print_CTU_Text(TComPic* pcPic, int intra_period)
{
#if CTU_INFO_CSV
    FILE* fp;
    char file_name[30];

    const TComSPS sps = pcPic->getPicSym()->getSPS();
    //const TComPPS pps = pcPic->getPicSym()->getPPS();
    const int frame_width = sps.getPicWidthInLumaSamples();
    const int frame_height = sps.getPicHeightInLumaSamples();
    const int frame_size = frame_width * frame_height;
    //const int ctu_size = sps.getMaxCUWidth();

    int x, sum_hdr, sum_data, sum_skip, sum_dist[3];
    const int frameno = encTopStat.m_encFrameNum;
    const int frame_sw = frameno & 0x1;
    const int frame_width_ctu = pcPic->getFrameWidthInCtus();
    const int frame_height_ctu = pcPic->getFrameHeightInCtus();
    const int ctu_num = pcPic->getNumberOfCtusInFrame();

    sprintf(file_name, "%s%d%s", "ctu_info/frame_", frameno, ".csv");
    fp = fopen(file_name, "w");

    sum_hdr = 0;
    sum_data = 0;
    sum_skip = 0;
    sum_dist[0] = 0;
    sum_dist[1] = 0;
    sum_dist[2] = 0;

    fprintf(fp, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n",
            "CtuNo", "(X Y)", "QP", "PSNRY", "PSNRU", "PSNRV", "LambdaIdx", "Lambda", "LumaW", "Intra", "Inter", "Skip", "H_Bits", "D_Bits", "A_Bits", "LeftBits",
            "HAD", "SAD", "DGR", "Non0", "Dist", "RDCost");

    for (x = 0; x < ctu_num; x++)
    {
        fprintf(fp, "%d, (%03d %03d), %d, %2.2f, %2.2f, %2.2f, %d, %f, %f, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
                x, x % frame_width_ctu, x / frame_width_ctu,
                encTopStat.CtuText[frame_sw][x].m_avgQP,
                encTopStat.CtuText[frame_sw][x].m_PSNR[0],
                encTopStat.CtuText[frame_sw][x].m_PSNR[1],
                encTopStat.CtuText[frame_sw][x].m_PSNR[2],
                encTopStat.CtuText[frame_sw][x].m_CtuLambdaIdx,
                encTopStat.CtuText[frame_sw][x].m_CtuLambda,
            #if CTU_WEIGHT == 1
                encTopStat.CtuText[frame_sw][x].hipa_weight,
            #elif CTU_WEIGHT == 2
                encTopStat.CtuText[frame_sw][x].lumo_weight,
            #elif CTU_WEIGHT == 3
                encTopStat.CtuText[frame_sw][x].var_weight,
            #endif
                encTopStat.CtuText[frame_sw][x].m_CtuType[2],
                encTopStat.CtuText[frame_sw][x].m_CtuType[1],
                encTopStat.CtuText[frame_sw][x].m_CtuType[0],
                encTopStat.CtuText[frame_sw][x].m_ctuHeaderBits,
                encTopStat.CtuText[frame_sw][x].m_ctuDataBits,
                encTopStat.CtuText[frame_sw][x].m_TotalBits,
                encTopStat.CtuText[frame_sw][x].m_CtuLeftBits,
                encTopStat.CtuText[frame_sw][x].m_ctuHAD,
                encTopStat.CtuText[frame_sw][x].m_ctuSAD,
                //encTopStat.CtuText[frame_sw][x].m_ctuDHAD,
                encTopStat.CtuText[frame_sw][x].m_ctuDGR,
                encTopStat.CtuText[frame_sw][x].m_NonZeroCoeffNum,
                encTopStat.CtuText[frame_sw][x].m_Dist[0],  // Y only
                encTopStat.CtuText[frame_sw][x].m_RdCost);
        sum_hdr += encTopStat.CtuText[frame_sw][x].m_ctuHeaderBits;
        sum_data += encTopStat.CtuText[frame_sw][x].m_ctuDataBits;
        sum_skip += encTopStat.CtuText[frame_sw][x].m_CtuType[0];
        sum_dist[0] += encTopStat.CtuText[frame_sw][x].m_Dist[0];
        sum_dist[1] += encTopStat.CtuText[frame_sw][x].m_Dist[1];
        sum_dist[2] += encTopStat.CtuText[frame_sw][x].m_Dist[2];
    }

    encTopStat.FrameStat[frameno].frameBitsHeader = sum_hdr;
    encTopStat.FrameStat[frameno].frameBitsData = sum_data;
    encTopStat.FrameStat[frameno].realSkip8x8Num = sum_skip;
    encTopStat.FrameStat[frameno].m_TotalDistortion[0] = (Int64)sum_dist[0] * 128 / frame_size;
    encTopStat.FrameStat[frameno].m_TotalDistortion[1] = (Int64)sum_dist[1] * 128 / frame_size;
    encTopStat.FrameStat[frameno].m_TotalDistortion[2] = (Int64)sum_dist[2] * 128 / frame_size;


#if 0   // Est Skip Block and Bits
    // estimate bits based on predicted skip blocks
    int is_intra = encTopStat.FrameStat[frameno].slice_type == 2;
    int intra_div = is_intra ? 5 : 1;   // ?
    int skip_bits = sum_skip / total_8x8;
    int nonSkip_bits = encTopStat.FrameStat[frameno].frameBits - skip_bits;
    int nonSkip_8x8s = total_frame_8x8 - sum_skip;
    if (nonSkip_8x8s == 0)  // all skip
        encTopStat.FrameStat[frameno].m_nonSkipBits8x8 = encTopStat.FrameStat[frameno - 1].m_nonSkipBits8x8;
    else
        encTopStat.FrameStat[frameno].m_nonSkipBits8x8 = encTopStat.FrameStat[frameno - 1].m_nonSkipBits8x8 + (float(nonSkip_bits) / nonSkip_8x8s / intra_div)/10;
        //encTopStat.FrameStat[frameno].m_nonSkipBits8x8 = float(nonSkip_bits) / nonSkip_8x8s / intra_div;

    if (is_intra)
        encTopStat.FrameStat[frameno].estframeBits = encTopStat.FrameStat[frameno].frameBits;
    else
    {
        float SAD_ratio = (encTopStat.FrameStat[frameno - 1].m_avgSAD==0) ? 1 : (float(encTopStat.FrameStat[frameno].m_avgSAD) / encTopStat.FrameStat[frameno-1].m_avgSAD);
        encTopStat.FrameStat[frameno].estframeBits = nonSkip_bits + int(sum_skip * SAD_ratio * encTopStat.FrameStat[frameno - 1].m_nonSkipBits8x8);
    }
#endif
    
    fprintf(fp, "\nTotal H_Bits=%d\n", sum_hdr);
    fprintf(fp, "Total D_Bits=%d\n", sum_data);
    fprintf(fp, "Total Sum Bits=%d\n", sum_hdr + sum_data);
    fprintf(fp, "Total Frame Bits=%d\n", encTopStat.FrameStat[frameno].frameBits);
    //fprintf(fp, "Est Frame Bits=%d\n", encTopStat.FrameStat[frameno].estframeBits);
    fprintf(fp, "Total NonZeroCoeffs_I=%d\n", encTopStat.FrameStat[frameno].m_TotalNonZeroNum_I);
    fprintf(fp, "Total NonZeroCoeffs_P=%d\n", encTopStat.FrameStat[frameno].m_TotalNonZeroNum_P);
    fprintf(fp, "Total NonZeroCoeffs_YUV=%d\n", encTopStat.FrameStat[frameno].m_TotalNonZeroNum_I + encTopStat.FrameStat[frameno].m_TotalNonZeroNum_P);
    fprintf(fp, "Total Distortion=%d\n", encTopStat.FrameStat[frameno].m_TotalDistortion[2]);
    fprintf(fp, "Total Skip=%d\n", sum_skip);
    
    fprintf(fp, "QPC=(%d %d) QPC_W=(%f %f)\n", encTopCRC.CRCQPC_delta[0], encTopCRC.CRCQPC_delta[1], encTopStat.FrameStat[frameno].blk_weight_c[0], encTopStat.FrameStat[frameno].blk_weight_c[1]);

    fprintf(fp, "Avg HAD=%d\n", encTopStat.FrameStat[frameno].m_avgHAD);
    fprintf(fp, "Avg SAD=%d\n", encTopStat.FrameStat[frameno].m_avgSAD);
    fprintf(fp, "Avg DGR=%d\n", encTopStat.FrameStat[frameno].m_avgDGR);
    fprintf(fp, "PSNR_Y=%f\n", encTopStat.FrameStat[frameno].m_encPSNR[0]);
    fprintf(fp, "PSNR_U=%f\n", encTopStat.FrameStat[frameno].m_encPSNR[1]);
    fprintf(fp, "PSNR_V=%f\n", encTopStat.FrameStat[frameno].m_encPSNR[2]);
    fprintf(fp, "PSNR_YUV=%f\n", encTopStat.FrameStat[frameno].m_encPSNR[3]);
    fprintf(fp, "MSSSIM=%f\n", encTopStat.FrameStat[frameno].m_encSSIM[3]);

    memset(encTopStat.CtuText[!frame_sw], 0, sizeof(EncCtuText) * ctu_num);
    // For intra frame, estimate from previous inter frame
    if ((encTopStat.m_encFrameNum % intra_period) == 0 && encTopStat.m_encFrameNum!=0)
    {
        encTopStat.CtuText[!frame_sw]->m_CtuType[0] = encTopStat.CtuText[frame_sw]->m_CtuType[0];
        encTopStat.CtuText[!frame_sw]->m_CtuType[1] = encTopStat.CtuText[frame_sw]->m_CtuType[1];
    }
    fclose(fp);

    printf("PIC=%d, SliceQP=%d, AvgQP=%d, RD_Lambda=%f, EstBits=%d, ActualBits=%d, CPB=%3.1f\n",
           frameno,
           encTopStat.FrameStat[frameno].frame_SliceQP,
           encTopStat.FrameStat[frameno].frame_AvgQP,
           encTopCRC.CRCPicLambda,
           encTopCRC.CRCPicTargetBits,
           encTopStat.FrameStat[frameno].frameBits,
           (float)encTopCRC.CRCCpbState * 100 / encTopCRC.CRCCpbSize
    );
#else
    int sum_hdr = 0;
    int sum_data = 0;
    int sum_skip = 0;
    int sum_dist[3] = { 0, 0, 0 };
    const int frameno = encTopStat.m_encFrameNum;
    const int frame_sw = frameno & 0x1;
    const int ctu_num = pcPic->getNumberOfCtusInFrame();
    const TComSPS sps = pcPic->getPicSym()->getSPS();
    const int frame_width = sps.getPicWidthInLumaSamples();
    const int frame_height = sps.getPicHeightInLumaSamples();
    const int frame_size = frame_width * frame_height;

    for (int x = 0; x < ctu_num; x++)
    {
        sum_hdr += encTopStat.CtuText[frame_sw][x].m_ctuHeaderBits;
        sum_data += encTopStat.CtuText[frame_sw][x].m_ctuDataBits;
        sum_skip += encTopStat.CtuText[frame_sw][x].m_CtuType[0];
        sum_dist[0] += encTopStat.CtuText[frame_sw][x].m_Dist[0];
        sum_dist[1] += encTopStat.CtuText[frame_sw][x].m_Dist[1];
        sum_dist[2] += encTopStat.CtuText[frame_sw][x].m_Dist[2];
    }

    encTopStat.FrameStat[frameno].frameBitsHeader = sum_hdr;
    encTopStat.FrameStat[frameno].frameBitsData = sum_data;
    encTopStat.FrameStat[frameno].realSkip8x8Num = sum_skip;
    encTopStat.FrameStat[frameno].m_TotalDistortion[0] = (Int64)sum_dist[0] * 128 / frame_size;
    encTopStat.FrameStat[frameno].m_TotalDistortion[1] = (Int64)sum_dist[1] * 128 / frame_size;
    encTopStat.FrameStat[frameno].m_TotalDistortion[2] = (Int64)sum_dist[2] * 128 / frame_size;

    printf("PIC=%d, SliceQP=%d, AvgQP=%d, RD_Lambda=%f, EstBits=%d, ActualBits=%d, DistY=%d, CPB=%3.1f\n",
           frameno,
           encTopStat.FrameStat[frameno].frame_SliceQP,
           encTopStat.FrameStat[frameno].frame_AvgQP,
           encTopCRC.CRCPicLambda,
           encTopCRC.CRCPicTargetBits,
           encTopStat.FrameStat[frameno].frameBits,
           encTopStat.FrameStat[frameno].m_TotalDistortion[0],
           (float)encTopCRC.CRCCpbState * 100 / encTopCRC.CRCCpbSize
    );
#endif
    memset(encTopStat.CtuText[0], 0, sizeof(EncCtuText)* ctu_num);
    memset(encTopStat.CtuText[1], 0, sizeof(EncCtuText)* ctu_num);
}


void Print_Frame_Stat()
{
    FILE* fp = fopen("enc_stat.csv", "w");
    int nb;

    if (fp == nullptr)
    {
        printf("ENC: open stat file failed!\n");
        return;
    }

    const int frame_num = encTopStat.m_encFrameNum;

    fprintf(fp, "\n%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n",
            "FrameNo", "POC", "Slice", "SliceQP", "AvgQP",
            "PSNR-YUV", "PSNR-PRV", "GopLeftBits", "EstBits", "T_Bits", "H_Bits", "D_Bits",
            "PicLambda", "HAD_a", "HAD_b",
            "AVG-HAD", "AVG-SAD", "AVG-DGR",
            "Non0_I", "Non0_P", "Non0_A", "SkipNum", "Distor", "MSE_A",
            "CpbState", "VBuf", "EstVBuf", "CpbOver"
    );

    for (nb = 0; nb < frame_num; nb++)
    {
        fprintf(fp, "%d, %d, %s, %d, %d, %f, %2.2f, %d, %d, %d, %d, %d, %f, %f, %f, %d, %d, %d, %d, %d, %d, %d, %d, %d, %3.1f, %d, %d, %d\n",
                nb,
                encTopStat.FrameStat[nb].frame_POC,
                (encTopStat.FrameStat[nb].slice_type == 2) ? "I" : (encTopStat.FrameStat[nb].slice_type == 1) ? "P" : "B",
                encTopStat.FrameStat[nb].frame_SliceQP,
                encTopStat.FrameStat[nb].frame_AvgQP,
                encTopStat.FrameStat[nb].m_encPSNR[3],
                (nb > 0) ? (encTopStat.FrameStat[nb].m_encPSNR[3] - encTopStat.FrameStat[nb-1].m_encPSNR[3]) : 0,
                //encTopStat.FrameStat[nb].m_encSSIM[3],
                encTopStat.FrameStat[nb].gop_left_bits,
                encTopStat.FrameStat[nb].target_bits,
                encTopStat.FrameStat[nb].frameBits,
                encTopStat.FrameStat[nb].frameBitsHeader,
                encTopStat.FrameStat[nb].frameBitsData,
                encTopStat.FrameStat[nb].pic_lambda,
                encTopStat.FrameStat[nb].pic_HAD_alpha,
                encTopStat.FrameStat[nb].pic_HAD_beta,
                encTopStat.FrameStat[nb].m_avgHAD,
                encTopStat.FrameStat[nb].m_avgSAD,
                encTopStat.FrameStat[nb].m_avgDGR,
                encTopStat.FrameStat[nb].m_TotalNonZeroNum_I,
                encTopStat.FrameStat[nb].m_TotalNonZeroNum_P,
                encTopStat.FrameStat[nb].m_TotalNonZeroNum_I + encTopStat.FrameStat[nb].m_TotalNonZeroNum_P,
                encTopStat.FrameStat[nb].realSkip8x8Num,
                encTopStat.FrameStat[nb].m_TotalDistortion[2],
                int(encTopStat.FrameStat[nb].m_encMSE[3]*128),
                (float)encTopStat.FrameStat[nb].cpb_state * 100 / encTopCRC.CRCCpbSize,
                encTopStat.FrameStat[nb].cpb_vbuf,
                encTopStat.FrameStat[nb].est_vbuf,
                encTopStat.FrameStat[nb].cpb_overflow
        );

        encTopStat.m_maxQP = max(encTopStat.FrameStat[nb].frame_maxQP, encTopStat.m_maxQP);
        encTopStat.m_minQP = min(encTopStat.FrameStat[nb].frame_minQP, encTopStat.m_minQP);
    }
    fclose(fp);
    printf("RC: Sequence Max QP = %d, MinQP = %d\n", encTopStat.m_maxQP, encTopStat.m_minQP);


    // VBR data
    FILE* fp_vbr = fopen("enc_vbr.csv", "w");

    if (fp_vbr == nullptr)
    {
        printf("ENC: open VBR stat file failed!\n");
        return;
    }
    fprintf(fp, "\n%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n",
            "FrameNo", "Slice", "SliceQP", "AvgQP", "PSNR-YUV", "PSNR-Y", "PSNR-U", "PSNR-V", 
            "GopLeftBits", "EstBits", "T_Bits", "AVG-SAD", "PicLambda", "Dist-YUV", "Dist-Y", "Distor-UV",
            "MSE_A", "MSE_Y", "MSE_U", "MSE_V", "VBuf"
    );

    for (nb = 0; nb < frame_num; nb++)
    {
        fprintf(fp, "%d, %s, %d, %d, %f, %f, %f, %f, %d, %d, %d, %d, %f, %d, %d, %d, %d, %d, %d, %d, %d\n",
                nb,
                (encTopStat.FrameStat[nb].slice_type == 2) ? "I" : (encTopStat.FrameStat[nb].slice_type == 1) ? "P" : "B",
                encTopStat.FrameStat[nb].frame_SliceQP,
                encTopStat.FrameStat[nb].frame_AvgQP,
                encTopStat.FrameStat[nb].m_encPSNR[3],
                encTopStat.FrameStat[nb].m_encPSNR[0],
                encTopStat.FrameStat[nb].m_encPSNR[1],
                encTopStat.FrameStat[nb].m_encPSNR[2],
                encTopStat.FrameStat[nb].gop_left_bits,
                encTopStat.FrameStat[nb].target_bits,
                encTopStat.FrameStat[nb].frameBits,
                encTopStat.FrameStat[nb].m_avgSAD,
                encTopStat.FrameStat[nb].pic_lambda,
                encTopStat.FrameStat[nb].m_TotalDistortion[2],
                encTopStat.FrameStat[nb].m_TotalDistortion[0],
                encTopStat.FrameStat[nb].m_TotalDistortion[1],
                int(encTopStat.FrameStat[nb].m_encMSE[3] * 128),
                int(encTopStat.FrameStat[nb].m_encMSE[0] * 128),
                int(encTopStat.FrameStat[nb].m_encMSE[1] * 128),
                int(encTopStat.FrameStat[nb].m_encMSE[2] * 128),
                encTopStat.FrameStat[nb].cpb_vbuf
        );
    }
    fclose(fp);
}