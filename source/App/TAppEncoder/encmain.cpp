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

/** \file     encmain.cpp
    \brief    Encoder application main
*/

#include <time.h>
#include <iostream>
#include "TAppEncTop.h"
#include "TAppCommon/program_options_lite.h"

//! \ingroup TAppEncoder
//! \{

#include "../Lib/TLibCommon/Debug.h"

// ====================================================================================================================
// Main function
// ====================================================================================================================

int main(int argc, char* argv[])
{
    double lambda[52];
    double mul[52];
    for (int QP = 0; QP <= 51; QP++)
    {
        // lambda = alpha * pow(2, (QP-12)/3), QP = 4.2005*log(lambda)+13.7122
        lambda[QP] = 0.57 * pow(2, (QP-12)/3);

        //lambda[QP] = int(lambda[QP] * 1024);
    }
    mul[51] = 65535 / lambda[51];
    mul[50] = 54662 / lambda[50];
    /*double a = log(2)*4.2005/3;
    double b = 13.7122 - 12*a;
    double c = 1-a;
    double alpha[51];
    double beta[51];
    double tmp[51];
    for (int QP = 0; QP <= 51; QP++)
    {
        // lambda = alpha * pow(2, (QP-12)/3), QP = 4.2005*log(lambda)+13.7122
        alpha[QP] = pow(2.71828, (c * QP - b) / 4.2005);

        double a0 = 0.61150518227251605;
        double b0 = 1.0070425145788775;
        double c0 = 1.0303737737998839;
        // lambda = a0 * pow(2, (c0*QP-12)/3)
        beta[QP] = a0 * pow(b0, QP);
        tmp[QP] = alpha[QP] / beta[QP];
    }*/
    /*{
        Double ratio = pow((51984 * 1.1) / 1516368, -0.5);
        Double CRCPicIAlpha = 11.249050 * ratio;
        int bitsFromP = int(2088960 * (pow(75.707094 / CRCPicIAlpha, 1 / -0.5)));

        Double ratio2 = pow((51984 * 1.2) / 1516368, -1.87);
        Double CRCPicIAlpha2 = 0.079863 * ratio2;
        int bitsFromP2 = int(2088960 * (pow(75.707094 / CRCPicIAlpha2, 1 / -1.87)));

        Double diffLambda0 = pow((Double)51984 / 57777, 1.0 * -0.5);
        Double diffLambda1 = pow((Double)51984 / 57777, 1.2 * -0.5);
        Double diffLambda2 = pow((Double)51984 / 57777, 1.5 * -0.5);

        Double old_beta = -0.5;
        Double old_alpha = 17.721224;
        int PicTargetBits = 57777;
        int PicActualBits = 5672;
        Double t1 = log((Double)PicTargetBits);
        Double t2 = log((Double)PicActualBits);

        Double bpp = 60000 / (Double)(1920*1088);
        Double Lambda = old_alpha * pow(bpp, old_beta);
        Double Lambda2 = old_alpha / bpp;
        int estQP = Int(4.2005 * log(Lambda) + 13.7122 + 0.5);
        int estQP2 = Int(4.2005 * log(Lambda2) + 13.7122 + 0.5);

        Double diffLambda = (old_beta) * (log((Double)PicTargetBits) - log((Double)PicActualBits));
        Double new_alpha = old_alpha * exp(diffLambda);
        new_alpha = new_alpha;
    }*/
  TAppEncTop  cTAppEncTop;

  // print information
  fprintf( stdout, "\n" );
  fprintf( stdout, "HM software: Encoder Version [%s] (including RExt)", NV_VERSION );
  fprintf( stdout, NVM_ONOS );
  fprintf( stdout, NVM_COMPILEDBY );
  fprintf( stdout, NVM_BITS );
  fprintf( stdout, "\n\n" );

  // create application encoder class
  cTAppEncTop.create();

  // parse configuration
  try
  {
    if(!cTAppEncTop.parseCfg( argc, argv ))
    {
      cTAppEncTop.destroy();
#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
      EnvVar::printEnvVar();
#endif
      return 1;
    }
  }
  catch (df::program_options_lite::ParseFailure &e)
  {
    std::cerr << "Error parsing option \""<< e.arg <<"\" with argument \""<< e.val <<"\"." << std::endl;
    return 1;
  }

#if PRINT_MACRO_VALUES
  printMacroSettings();
#endif

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
  EnvVar::printEnvVarInUse();
#endif

  // starting time
  Double dResult;
  clock_t lBefore = clock();

  // call encoding function
  cTAppEncTop.encode();

  // ending time
  dResult = (Double)(clock()-lBefore) / CLOCKS_PER_SEC;
  printf("\n Total Time: %12.3f sec.\n", dResult);

  // destroy application encoder class
  cTAppEncTop.destroy();

  return 0;
}

//! \}
