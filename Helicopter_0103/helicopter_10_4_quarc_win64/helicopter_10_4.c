/*
 * helicopter_10_4.c
 *
 * Code generation for model "helicopter_10_4".
 *
 * Model version              : 1.168
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Wed Feb 22 17:32:35 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter_10_4.h"
#include "helicopter_10_4_private.h"
#include "helicopter_10_4_dt.h"

/* Block signals (auto storage) */
B_helicopter_10_4_T helicopter_10_4_B;

/* Continuous states */
X_helicopter_10_4_T helicopter_10_4_X;

/* Block states (auto storage) */
DW_helicopter_10_4_T helicopter_10_4_DW;

/* Real-time model */
RT_MODEL_helicopter_10_4_T helicopter_10_4_M_;
RT_MODEL_helicopter_10_4_T *const helicopter_10_4_M = &helicopter_10_4_M_;

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helicopter_10_4_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter_10_4_output(void)
{
  /* local block i/o variables */
  real_T rtb_Sum3_o[2];
  real_T rtb_Sum2_f[6];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T rtb_Derivative;
  int32_T i;
  int32_T i_0;
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
    /* set solver stop time */
    if (!(helicopter_10_4_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter_10_4_M->solverInfo,
                            ((helicopter_10_4_M->Timing.clockTickH0 + 1) *
        helicopter_10_4_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter_10_4_M->solverInfo,
                            ((helicopter_10_4_M->Timing.clockTick0 + 1) *
        helicopter_10_4_M->Timing.stepSize0 +
        helicopter_10_4_M->Timing.clockTickH0 *
        helicopter_10_4_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter_10_4_M)) {
    helicopter_10_4_M->Timing.t[0] = rtsiGetT(&helicopter_10_4_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter_10_4/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicopter_10_4_DW.HILReadEncoderTimebase_Task, 1,
         &helicopter_10_4_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter_10_4_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter_10_4_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter_10_4_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S4>/Travel: Count to rad' */
    helicopter_10_4_B.TravelCounttorad = helicopter_10_4_P.TravelCounttorad_Gain
      * rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S11>/Gain' */
    helicopter_10_4_B.Gain = helicopter_10_4_P.Gain_Gain *
      helicopter_10_4_B.TravelCounttorad;

    /* Sum: '<Root>/Sum1' incorporates:
     *  Constant: '<Root>/ Travel_offset'
     */
    helicopter_10_4_B.Sum1 = helicopter_10_4_P.Travel_offset_Value +
      helicopter_10_4_B.Gain;
  }

  /* TransferFcn: '<S4>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter_10_4_P.TravelTransferFcn_C *
    helicopter_10_4_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helicopter_10_4_P.TravelTransferFcn_D *
    helicopter_10_4_B.TravelCounttorad;

  /* Gain: '<S12>/Gain' */
  helicopter_10_4_B.Gain_d = helicopter_10_4_P.Gain_Gain_l * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter_10_4_B.PitchCounttorad = helicopter_10_4_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S8>/Gain' */
    helicopter_10_4_B.Gain_i = helicopter_10_4_P.Gain_Gain_a *
      helicopter_10_4_B.PitchCounttorad;
  }

  /* TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter_10_4_P.PitchTransferFcn_C *
    helicopter_10_4_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helicopter_10_4_P.PitchTransferFcn_D *
    helicopter_10_4_B.PitchCounttorad;

  /* Gain: '<S9>/Gain' */
  helicopter_10_4_B.Gain_b = helicopter_10_4_P.Gain_Gain_ae * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    helicopter_10_4_B.ElevationCounttorad =
      helicopter_10_4_P.ElevationCounttorad_Gain * rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S6>/Gain' */
    helicopter_10_4_B.Gain_e = helicopter_10_4_P.Gain_Gain_lv *
      helicopter_10_4_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter_10_4_B.Sum = helicopter_10_4_B.Gain_e +
      helicopter_10_4_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter_10_4_P.ElevationTransferFcn_C *
    helicopter_10_4_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helicopter_10_4_P.ElevationTransferFcn_D *
    helicopter_10_4_B.ElevationCounttorad;

  /* Gain: '<S7>/Gain' */
  helicopter_10_4_B.Gain_dg = helicopter_10_4_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  helicopter_10_4_B.Gain1[0] = helicopter_10_4_P.Gain1_Gain *
    helicopter_10_4_B.Sum1;
  helicopter_10_4_B.Gain1[1] = helicopter_10_4_P.Gain1_Gain *
    helicopter_10_4_B.Gain_d;
  helicopter_10_4_B.Gain1[2] = helicopter_10_4_P.Gain1_Gain *
    helicopter_10_4_B.Gain_i;
  helicopter_10_4_B.Gain1[3] = helicopter_10_4_P.Gain1_Gain *
    helicopter_10_4_B.Gain_b;
  helicopter_10_4_B.Gain1[4] = helicopter_10_4_P.Gain1_Gain *
    helicopter_10_4_B.Sum;
  helicopter_10_4_B.Gain1[5] = helicopter_10_4_P.Gain1_Gain *
    helicopter_10_4_B.Gain_dg;
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter_10_4_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter_10_4_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_10_4_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter_10_4_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[104]) {
      currTimeIndex = 103;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter_10_4_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_Sum3_o[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 105;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_Sum3_o[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 105;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 2; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum3_o[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 105;
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helicopter_10_4_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter_10_4_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_10_4_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helicopter_10_4_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[104]) {
      currTimeIndex = 103;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter_10_4_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum2_f[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 105;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum2_f[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 105;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 6; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum2_f[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 105;
          }
        }
      }
    }
  }

  /* Sum: '<Root>/Sum2' */
  for (i = 0; i < 6; i++) {
    rtb_Sum2_f[i] = helicopter_10_4_B.Gain1[i] - rtb_Sum2_f[i];
  }

  /* End of Sum: '<Root>/Sum2' */

  /* Sum: '<Root>/Sum3' incorporates:
   *  Gain: '<Root>/Gain'
   */
  for (i = 0; i < 2; i++) {
    rtb_Derivative = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      rtb_Derivative += helicopter_10_4_P.K[(i_0 << 1) + i] * rtb_Sum2_f[i_0];
    }

    rtb_Sum3_o[i] -= rtb_Derivative;
  }

  /* End of Sum: '<Root>/Sum3' */

  /* Sum: '<S5>/Sum' incorporates:
   *  Constant: '<S5>/Vd_bias'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<S5>/Sum2'
   *  Sum: '<S5>/Sum3'
   */
  helicopter_10_4_B.Sum_k = ((rtb_Sum3_o[0] - helicopter_10_4_B.Gain1[2]) *
    helicopter_10_4_P.K_pp - helicopter_10_4_P.K_pd * helicopter_10_4_B.Gain1[3])
    + helicopter_10_4_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
  }

  /* Integrator: '<S3>/Integrator'
   *
   * Regarding '<S3>/Integrator':
   *  Limited Integrator
   */
  if (helicopter_10_4_X.Integrator_CSTATE >=
      helicopter_10_4_P.Integrator_UpperSat ) {
    helicopter_10_4_X.Integrator_CSTATE = helicopter_10_4_P.Integrator_UpperSat;
  } else if (helicopter_10_4_X.Integrator_CSTATE <=
             (helicopter_10_4_P.Integrator_LowerSat) ) {
    helicopter_10_4_X.Integrator_CSTATE = (helicopter_10_4_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter_10_4_X.Integrator_CSTATE;

  /* Sum: '<S3>/Sum' */
  rtb_Derivative = rtb_Sum3_o[1] - helicopter_10_4_B.Gain1[4];

  /* Sum: '<S3>/Sum2' incorporates:
   *  Constant: '<S3>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  helicopter_10_4_B.Sum2 = ((helicopter_10_4_P.K_ep * rtb_Derivative +
    rtb_Backgain) - helicopter_10_4_P.K_ed * helicopter_10_4_B.Gain1[5]) +
    helicopter_10_4_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helicopter_10_4_B.Sum2 - helicopter_10_4_B.Sum_k) *
    helicopter_10_4_P.Backgain_Gain;

  /* Gain: '<S3>/K_ei' */
  helicopter_10_4_B.K_ei = helicopter_10_4_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helicopter_10_4_DW.TimeStampA >= helicopter_10_4_M->Timing.t[0]) &&
      (helicopter_10_4_DW.TimeStampB >= helicopter_10_4_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Derivative = helicopter_10_4_DW.TimeStampA;
    lastU = &helicopter_10_4_DW.LastUAtTimeA;
    if (helicopter_10_4_DW.TimeStampA < helicopter_10_4_DW.TimeStampB) {
      if (helicopter_10_4_DW.TimeStampB < helicopter_10_4_M->Timing.t[0]) {
        rtb_Derivative = helicopter_10_4_DW.TimeStampB;
        lastU = &helicopter_10_4_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter_10_4_DW.TimeStampA >= helicopter_10_4_M->Timing.t[0]) {
        rtb_Derivative = helicopter_10_4_DW.TimeStampB;
        lastU = &helicopter_10_4_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter_10_4_B.PitchCounttorad - *lastU) /
      (helicopter_10_4_M->Timing.t[0] - rtb_Derivative);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S10>/Gain' */
  helicopter_10_4_B.Gain_l = helicopter_10_4_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
  }

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter_10_4_P.BackmotorSaturation_UpperSat) {
    helicopter_10_4_B.BackmotorSaturation =
      helicopter_10_4_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helicopter_10_4_P.BackmotorSaturation_LowerSat) {
    helicopter_10_4_B.BackmotorSaturation =
      helicopter_10_4_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter_10_4_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Derivative = (helicopter_10_4_B.Sum_k + helicopter_10_4_B.Sum2) *
    helicopter_10_4_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Derivative > helicopter_10_4_P.FrontmotorSaturation_UpperSat) {
    helicopter_10_4_B.FrontmotorSaturation =
      helicopter_10_4_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Derivative < helicopter_10_4_P.FrontmotorSaturation_LowerSat) {
    helicopter_10_4_B.FrontmotorSaturation =
      helicopter_10_4_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter_10_4_B.FrontmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter_10_4/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter_10_4_DW.HILWriteAnalog_Buffer[0] =
        helicopter_10_4_B.FrontmotorSaturation;
      helicopter_10_4_DW.HILWriteAnalog_Buffer[1] =
        helicopter_10_4_B.BackmotorSaturation;
      result = hil_write_analog(helicopter_10_4_DW.HILInitialize_Card,
        helicopter_10_4_P.HILWriteAnalog_channels, 2,
        &helicopter_10_4_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter_10_4_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter_10_4_DW.TimeStampA == (rtInf)) {
    helicopter_10_4_DW.TimeStampA = helicopter_10_4_M->Timing.t[0];
    lastU = &helicopter_10_4_DW.LastUAtTimeA;
  } else if (helicopter_10_4_DW.TimeStampB == (rtInf)) {
    helicopter_10_4_DW.TimeStampB = helicopter_10_4_M->Timing.t[0];
    lastU = &helicopter_10_4_DW.LastUAtTimeB;
  } else if (helicopter_10_4_DW.TimeStampA < helicopter_10_4_DW.TimeStampB) {
    helicopter_10_4_DW.TimeStampA = helicopter_10_4_M->Timing.t[0];
    lastU = &helicopter_10_4_DW.LastUAtTimeA;
  } else {
    helicopter_10_4_DW.TimeStampB = helicopter_10_4_M->Timing.t[0];
    lastU = &helicopter_10_4_DW.LastUAtTimeB;
  }

  *lastU = helicopter_10_4_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter_10_4_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter_10_4_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicopter_10_4_M->Timing.clockTick0)) {
    ++helicopter_10_4_M->Timing.clockTickH0;
  }

  helicopter_10_4_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter_10_4_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helicopter_10_4_M->Timing.clockTick1)) {
      ++helicopter_10_4_M->Timing.clockTickH1;
    }

    helicopter_10_4_M->Timing.t[1] = helicopter_10_4_M->Timing.clockTick1 *
      helicopter_10_4_M->Timing.stepSize1 +
      helicopter_10_4_M->Timing.clockTickH1 *
      helicopter_10_4_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter_10_4_derivatives(void)
{
  XDot_helicopter_10_4_T *_rtXdot;
  _rtXdot = ((XDot_helicopter_10_4_T *) helicopter_10_4_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_10_4_P.TravelTransferFcn_A *
    helicopter_10_4_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_10_4_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_10_4_P.PitchTransferFcn_A *
    helicopter_10_4_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_10_4_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter_10_4_P.ElevationTransferFcn_A *
    helicopter_10_4_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter_10_4_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter_10_4_X.Integrator_CSTATE <=
            (helicopter_10_4_P.Integrator_LowerSat) );
    usat = ( helicopter_10_4_X.Integrator_CSTATE >=
            helicopter_10_4_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter_10_4_B.K_ei > 0)) ||
        (usat && (helicopter_10_4_B.K_ei < 0)) ) {
      ((XDot_helicopter_10_4_T *) helicopter_10_4_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter_10_4_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter_10_4_T *) helicopter_10_4_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter_10_4_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_10_4/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter_10_4_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helicopter_10_4_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter_10_4_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
      return;
    }

    if ((helicopter_10_4_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter_10_4_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter_10_4_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = helicopter_10_4_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter_10_4_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicopter_10_4_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(helicopter_10_4_DW.HILInitialize_Card,
        helicopter_10_4_P.HILInitialize_analog_input_chan, 8U,
        &helicopter_10_4_DW.HILInitialize_AIMinimums[0],
        &helicopter_10_4_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_10_4_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter_10_4_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter_10_4_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = helicopter_10_4_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter_10_4_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicopter_10_4_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges
        (helicopter_10_4_DW.HILInitialize_Card,
         helicopter_10_4_P.HILInitialize_analog_output_cha, 8U,
         &helicopter_10_4_DW.HILInitialize_AOMinimums[0],
         &helicopter_10_4_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_10_4_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter_10_4_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_10_4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_10_4_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter_10_4_DW.HILInitialize_Card,
        helicopter_10_4_P.HILInitialize_analog_output_cha, 8U,
        &helicopter_10_4_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_10_4_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_10_4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_10_4_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter_10_4_DW.HILInitialize_Card,
         helicopter_10_4_P.HILInitialize_analog_output_cha, 8U,
         &helicopter_10_4_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_10_4_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter_10_4_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter_10_4_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter_10_4_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter_10_4_DW.HILInitialize_Card,
         helicopter_10_4_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter_10_4_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_10_4_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter_10_4_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter_10_4_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter_10_4_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter_10_4_DW.HILInitialize_Card,
        helicopter_10_4_P.HILInitialize_encoder_channels, 8U,
        &helicopter_10_4_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_10_4_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter_10_4_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_10_4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_10_4_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter_10_4_DW.HILInitialize_Card,
        helicopter_10_4_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter_10_4_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter_10_4_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues =
          &helicopter_10_4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter_10_4_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter_10_4_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = helicopter_10_4_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter_10_4_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = p_HILInitialize_pwm_channels[i1];
            helicopter_10_4_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              helicopter_10_4_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter_10_4_DW.HILInitialize_Card,
          &helicopter_10_4_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes, &helicopter_10_4_DW.HILInitialize_POSortedFreqs
          [0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter_10_4_DW.HILInitialize_Card,
          &helicopter_10_4_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter_10_4_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_10_4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            helicopter_10_4_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter_10_4_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter_10_4_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter_10_4_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter_10_4_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter_10_4_DW.HILInitialize_Card,
        helicopter_10_4_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &helicopter_10_4_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter_10_4_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter_10_4_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicopter_10_4_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] =
            helicopter_10_4_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter_10_4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_10_4_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter_10_4_DW.HILInitialize_Card,
        helicopter_10_4_P.HILInitialize_pwm_channels, 8U,
        &helicopter_10_4_DW.HILInitialize_POSortedFreqs[0],
        &helicopter_10_4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_10_4_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter_10_4_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_10_4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_10_4_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter_10_4_DW.HILInitialize_Card,
        helicopter_10_4_P.HILInitialize_pwm_channels, 8U,
        &helicopter_10_4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_10_4_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_10_4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_10_4_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter_10_4_DW.HILInitialize_Card,
         helicopter_10_4_P.HILInitialize_pwm_channels, 8U,
         &helicopter_10_4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter_10_4/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (helicopter_10_4_DW.HILInitialize_Card,
       helicopter_10_4_P.HILReadEncoderTimebase_samples_,
       helicopter_10_4_P.HILReadEncoderTimebase_channels, 3,
       &helicopter_10_4_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -8.1603049590754432E-8,
      -3.7650934206307147E-7, -6.2155353888674982E-7, -4.8848677632778625E-7,
      -2.7670414189540363E-7, -2.9658893489262082E-7, -5.28595818855844E-7,
      -7.90899827572553E-7, -1.2056622605312566E-6, -1.5722765146058694E-6,
      -1.7315481394276065E-6, -1.5925362263635215E-6, -1.4231301547614639E-6,
      -1.3421743199858339E-6, -1.3160829247676137E-6, -1.1109343347480953E-6,
      -6.954085262742544E-7, -1.1578597577429036E-7, 2.6528675727558614E-7,
      1.5853359338049751E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.00025401352066387137, 0.00045630215641509952, 0.000467967276828556,
      0.00045686181902391773, 0.00045844028952447125, 0.00045858338526347865,
      0.00045850118554257755, 0.0004585088363017112, 0.00045851208872685155,
      0.00045851290657504106, 0.00045851381038110431, 0.00045851741880156621,
      0.00045850803787707251, 0.00045851217578095859, 0.00045885525324053454,
      0.00045729723260323107, 0.00045563288756372285, 0.00050696897159787529,
      0.00030463130032885231, 6.3509449328405868E-5, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

    helicopter_10_4_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_10_4_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_10_4_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace1' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, -3.1415926535892149, -2.0794002176192797E-26,
      3.9472984687604958E-23, 3.7409191866524359E-8, 3.7929652528632605E-8,
      2.71855409040797E-8, 3.1764077859042492E-9, 2.7332939936626098E-8,
      2.7274778206639405E-8, 5.5575644381263E-8, 7.412869421529721E-8,
      1.0578197197622695E-7, 9.9784840053575514E-8, 9.1832602221814541E-8,
      7.6164294937913786E-8, 8.0674821798868231E-8, 7.80897507736625E-8,
      6.0820851230493548E-8, 3.7318488294813007E-8, -1.2610787542677162E-8,
      -2.2433435646550535E-8, 4.2179281141895472E-12, -7.2253364305147541E-12,
      1.2345124058010788E-11, -2.1032356574104446E-11, 3.5718281717759665E-11,
      -6.0440978884234448E-11, 1.0186029121801326E-10, -1.7086827679986133E-10,
      2.8509845293301562E-10, -4.727412137272357E-10, 7.7814938191519461E-10,
      -1.2696562424183956E-9, 2.04955847530102E-9, -3.2647242856481062E-9,
      5.1124379683577256E-9, -7.8270700005877386E-9, 1.1613152984637708E-8,
      -1.6447643942790291E-8, 2.1584902523819858E-8, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, -1.3916693562943657E-22, 3.7409191866523538E-8, 3.7929652528631996E-8,
      2.7185540904079705E-8, 3.176407785904044E-9, 2.7332939936626204E-8,
      2.7274778206639259E-8, 5.5575644381263135E-8, 7.4128694215297475E-8,
      1.0578197197622489E-7, 9.9784840053573979E-8, 9.1832602221814488E-8,
      7.6164294937914475E-8, 8.0674821798868178E-8, 7.8089750773662609E-8,
      6.0820851230492966E-8, 3.7318488294813643E-8, -1.2610787542677421E-8,
      -2.2433435646550747E-8, 4.2179281142820977E-12, -7.2253364305143276E-12,
      1.2345124058005156E-11, -2.1032356573993623E-11, 3.5718281717747968E-11,
      -6.0440978884256084E-11, 1.0186029121801529E-10, -1.708682767998627E-10,
      2.8509845293301681E-10, -4.7274121372721005E-10, 7.7814938191519461E-10,
      -1.2696562424183919E-9, 2.0495584753010214E-9, -3.2647242856481095E-9,
      5.1124379683577264E-9, -7.82707000058789E-9, 1.1613152984637518E-8,
      -1.6447643942790222E-8, 2.1584902523819458E-8, -2.4415282711109631E-8, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 3.0580411726565181E-22, -6.6098470168508892E-8,
      -6.7018074464453763E-8, -4.8034255080263086E-8, -5.6124092717303191E-9,
      -4.8294695096995677E-8, -4.8191928873436965E-8, -9.8196857214600615E-8,
      -1.3097832481126511E-7, -1.8690664433447354E-7, -1.7631028483812231E-7,
      -1.6225943987544757E-7, -1.3457503692734112E-7, -1.4254470722191891E-7,
      -1.3797713354503939E-7, -1.0746463690063528E-7, -6.5938205617747131E-8,
      2.2282057499804002E-8, 3.9637738825030131E-8, -7.4526762465458888E-12,
      1.2766479591871067E-11, -2.1812655488225179E-11, 3.7162165880039451E-11,
      -6.3110793385056376E-11, 1.0679343873543145E-10, -1.7997740887995144E-10,
      3.0190793046526048E-10, -5.037417449037357E-10, 8.3528858694590633E-10,
      -1.3749156595171554E-9, 2.2433613525571689E-9, -3.6213741323701474E-9,
      5.7684560942474375E-9, -9.0331897504112E-9, 1.3829685356118602E-8,
      -2.05193324906949E-8, 2.9061416412677249E-8, -3.8138461815786933E-8,
      4.3139473359758751E-8, -3.1733487811982211E-8, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -6.60984701685108E-8, -6.7018074464453234E-8, -4.8034255080263026E-8,
      -5.6124092717304051E-9, -4.8294695096995174E-8, -4.8191928873437018E-8,
      -9.8196857214601409E-8, -1.3097832481126119E-7, -1.8690664433447915E-7,
      -1.7631028483811797E-7, -1.6225943987544609E-7, -1.3457503692733974E-7,
      -1.4254470722192208E-7, -1.3797713354503505E-7, -1.0746463690063607E-7,
      -6.5938205617747713E-8, 2.2282057499804095E-8, 3.9637738825029364E-8,
      -7.4526762465975359E-12, 1.2766479591356375E-11, -2.1812655488152769E-11,
      3.7162165880066838E-11, -6.311079338493416E-11, 1.067934387355165E-10,
      -1.7997740887994976E-10, 3.0190793046526803E-10, -5.03741744903724E-10,
      8.3528858694588968E-10, -1.3749156595172058E-9, 2.2433613525571912E-9,
      -3.6213741323701941E-9, 5.7684560942474292E-9, -9.0331897504112153E-9,
      1.3829685356118602E-8, -2.0519332490694915E-8, 2.9061416412677249E-8,
      -3.8138461815786986E-8, 4.31394733597588E-8, -3.1733487811982178E-8,
      -2.5531337562482489E-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.5918191725561824E-28,
      1.5875845041491961E-5, 1.2643039734451764E-5, 1.2635953806960002E-5,
      1.2757149948421913E-5, 1.2736379695117546E-5, 1.2735794396744398E-5,
      1.2736434775887321E-5, 1.2736418893783526E-5, 1.2736477957672868E-5,
      1.2736473979821318E-5, 1.2736519679579481E-5, 1.2736700500563003E-5,
      1.2735921946859126E-5, 1.2736913914310041E-5, 1.2757558926508564E-5,
      1.2639289632615858E-5, 1.2648376108489708E-5, 1.5877362208223569E-5,
      3.5207276843611317E-14, -4.2378021268493794E-15, -4.5640170840513822E-15,
      5.6234676157604488E-15, -4.4824633447501031E-15, 3.0765964408088268E-15,
      -1.9559806046251926E-15, 1.1868314944223908E-15, -6.9783634326699159E-16,
      4.011284696631098E-16, -2.266693838473415E-16, 1.2638726643247762E-16,
      -6.9719920470645785E-17, 3.8123103864502963E-17, -2.0693123747097604E-17,
      1.1162347780835729E-17, -5.9890668432449379E-18, 3.1984798988110419E-18,
      -1.7012131876211195E-18, 9.015932153453259E-19, -4.7628991482143181E-19,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.5875845041491957E-5, 1.2643039734451761E-5,
      1.2635953806960005E-5, 1.2757149948421916E-5, 1.2736379695117541E-5,
      1.2735794396744398E-5, 1.273643477588732E-5, 1.273641889378353E-5,
      1.2736477957672866E-5, 1.2736473979821322E-5, 1.2736519679579481E-5,
      1.2736700500563086E-5, 1.2735921946859077E-5, 1.2736913914310081E-5,
      1.2757558926508614E-5, 1.2639289632615812E-5, 1.2648376108489735E-5,
      1.5877362208223573E-5, 3.520727684361458E-14, -4.237802126850207E-15,
      -4.5640170840503567E-15, 5.6234676157614656E-15, -4.4824633447488228E-15,
      3.0765964408105374E-15, -1.9559806046251658E-15, 1.186831494423188E-15,
      -6.9783634326752861E-16, 4.0112846966421431E-16, -2.2666938384804565E-16,
      1.2638726643200665E-16, -6.9719920471538E-17, 3.8123103863653027E-17,
      -2.0693123746689809E-17, 1.1162347780582646E-17, -5.989066844232779E-18,
      3.198479898802228E-18, -1.7012131885194866E-18, 9.01593214049836E-19,
      -4.7628991634109453E-19, 2.5089612580797149E-19, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter_10_4_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_10_4_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_10_4_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter_10_4_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter_10_4_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter_10_4_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter_10_4_X.Integrator_CSTATE = helicopter_10_4_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter_10_4_DW.TimeStampA = (rtInf);
  helicopter_10_4_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter_10_4_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_10_4/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter_10_4_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter_10_4_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter_10_4_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter_10_4_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_10_4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_10_4_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter_10_4_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter_10_4_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_10_4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_10_4_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter_10_4_DW.HILInitialize_Card
                         , helicopter_10_4_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter_10_4_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter_10_4_DW.HILInitialize_AOVoltages[0]
                         , &helicopter_10_4_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicopter_10_4_DW.HILInitialize_Card,
            helicopter_10_4_P.HILInitialize_analog_output_cha,
            num_final_analog_outputs,
            &helicopter_10_4_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter_10_4_DW.HILInitialize_Card,
            helicopter_10_4_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &helicopter_10_4_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_10_4_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter_10_4_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter_10_4_DW.HILInitialize_Card);
    hil_close(helicopter_10_4_DW.HILInitialize_Card);
    helicopter_10_4_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  helicopter_10_4_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter_10_4_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  helicopter_10_4_initialize();
}

void MdlTerminate(void)
{
  helicopter_10_4_terminate();
}

/* Registration function */
RT_MODEL_helicopter_10_4_T *helicopter_10_4(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter_10_4_P.Integrator_UpperSat = rtInf;
  helicopter_10_4_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter_10_4_M, 0,
                sizeof(RT_MODEL_helicopter_10_4_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter_10_4_M->solverInfo,
                          &helicopter_10_4_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter_10_4_M->solverInfo, &rtmGetTPtr(helicopter_10_4_M));
    rtsiSetStepSizePtr(&helicopter_10_4_M->solverInfo,
                       &helicopter_10_4_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter_10_4_M->solverInfo,
                 &helicopter_10_4_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter_10_4_M->solverInfo, (real_T **)
                         &helicopter_10_4_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter_10_4_M->solverInfo,
      &helicopter_10_4_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter_10_4_M->solverInfo, (&rtmGetErrorStatus
      (helicopter_10_4_M)));
    rtsiSetRTModelPtr(&helicopter_10_4_M->solverInfo, helicopter_10_4_M);
  }

  rtsiSetSimTimeStep(&helicopter_10_4_M->solverInfo, MAJOR_TIME_STEP);
  helicopter_10_4_M->ModelData.intgData.f[0] = helicopter_10_4_M->
    ModelData.odeF[0];
  helicopter_10_4_M->ModelData.contStates = ((real_T *) &helicopter_10_4_X);
  rtsiSetSolverData(&helicopter_10_4_M->solverInfo, (void *)
                    &helicopter_10_4_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter_10_4_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter_10_4_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter_10_4_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter_10_4_M->Timing.sampleTimes =
      (&helicopter_10_4_M->Timing.sampleTimesArray[0]);
    helicopter_10_4_M->Timing.offsetTimes =
      (&helicopter_10_4_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter_10_4_M->Timing.sampleTimes[0] = (0.0);
    helicopter_10_4_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter_10_4_M->Timing.offsetTimes[0] = (0.0);
    helicopter_10_4_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter_10_4_M, &helicopter_10_4_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter_10_4_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter_10_4_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter_10_4_M, 25.0);
  helicopter_10_4_M->Timing.stepSize0 = 0.002;
  helicopter_10_4_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter_10_4_M->Sizes.checksums[0] = (4056607422U);
  helicopter_10_4_M->Sizes.checksums[1] = (2765574339U);
  helicopter_10_4_M->Sizes.checksums[2] = (795282245U);
  helicopter_10_4_M->Sizes.checksums[3] = (3561673981U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter_10_4_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter_10_4_M->extModeInfo,
      &helicopter_10_4_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter_10_4_M->extModeInfo,
                        helicopter_10_4_M->Sizes.checksums);
    rteiSetTPtr(helicopter_10_4_M->extModeInfo, rtmGetTPtr(helicopter_10_4_M));
  }

  helicopter_10_4_M->solverInfoPtr = (&helicopter_10_4_M->solverInfo);
  helicopter_10_4_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter_10_4_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter_10_4_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter_10_4_M->ModelData.blockIO = ((void *) &helicopter_10_4_B);

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      helicopter_10_4_B.Gain1[i] = 0.0;
    }

    helicopter_10_4_B.TravelCounttorad = 0.0;
    helicopter_10_4_B.Gain = 0.0;
    helicopter_10_4_B.Sum1 = 0.0;
    helicopter_10_4_B.Gain_d = 0.0;
    helicopter_10_4_B.PitchCounttorad = 0.0;
    helicopter_10_4_B.Gain_i = 0.0;
    helicopter_10_4_B.Gain_b = 0.0;
    helicopter_10_4_B.ElevationCounttorad = 0.0;
    helicopter_10_4_B.Gain_e = 0.0;
    helicopter_10_4_B.Sum = 0.0;
    helicopter_10_4_B.Gain_dg = 0.0;
    helicopter_10_4_B.Sum_k = 0.0;
    helicopter_10_4_B.Sum2 = 0.0;
    helicopter_10_4_B.K_ei = 0.0;
    helicopter_10_4_B.Gain_l = 0.0;
    helicopter_10_4_B.BackmotorSaturation = 0.0;
    helicopter_10_4_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter_10_4_M->ModelData.defaultParam = ((real_T *)&helicopter_10_4_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter_10_4_X;
    helicopter_10_4_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter_10_4_X, 0,
                  sizeof(X_helicopter_10_4_T));
  }

  /* states (dwork) */
  helicopter_10_4_M->ModelData.dwork = ((void *) &helicopter_10_4_DW);
  (void) memset((void *)&helicopter_10_4_DW, 0,
                sizeof(DW_helicopter_10_4_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_10_4_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_10_4_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_10_4_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_10_4_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_10_4_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_10_4_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_10_4_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_10_4_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter_10_4_DW.TimeStampA = 0.0;
  helicopter_10_4_DW.LastUAtTimeA = 0.0;
  helicopter_10_4_DW.TimeStampB = 0.0;
  helicopter_10_4_DW.LastUAtTimeB = 0.0;
  helicopter_10_4_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter_10_4_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter_10_4_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter_10_4_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter_10_4_M->Sizes.numY = (0); /* Number of model outputs */
  helicopter_10_4_M->Sizes.numU = (0); /* Number of model inputs */
  helicopter_10_4_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter_10_4_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter_10_4_M->Sizes.numBlocks = (59);/* Number of blocks */
  helicopter_10_4_M->Sizes.numBlockIO = (18);/* Number of block outputs */
  helicopter_10_4_M->Sizes.numBlockPrms = (153);/* Sum of parameter "widths" */
  return helicopter_10_4_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
