/*
 * helicopter_statefeedback.c
 *
 * Code generation for model "helicopter_statefeedback".
 *
 * Model version              : 1.166
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Wed Feb 22 13:47:52 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter_statefeedback.h"
#include "helicopter_statefeedback_private.h"
#include "helicopter_statefeedback_dt.h"

/* Block signals (auto storage) */
B_helicopter_statefeedback_T helicopter_statefeedback_B;

/* Continuous states */
X_helicopter_statefeedback_T helicopter_statefeedback_X;

/* Block states (auto storage) */
DW_helicopter_statefeedback_T helicopter_statefeedback_DW;

/* Real-time model */
RT_MODEL_helicopter_statefeed_T helicopter_statefeedback_M_;
RT_MODEL_helicopter_statefeed_T *const helicopter_statefeedback_M =
  &helicopter_statefeedback_M_;

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
  helicopter_statefeedback_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter_statefeedback_output(void)
{
  /* local block i/o variables */
  real_T rtb_Sum2_f[4];
  real_T rtb_Backgain;
  real_T rtb_Derivative;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T lastTime;
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
    /* set solver stop time */
    if (!(helicopter_statefeedback_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter_statefeedback_M->solverInfo,
                            ((helicopter_statefeedback_M->Timing.clockTickH0 + 1)
        * helicopter_statefeedback_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter_statefeedback_M->solverInfo,
                            ((helicopter_statefeedback_M->Timing.clockTick0 + 1)
        * helicopter_statefeedback_M->Timing.stepSize0 +
        helicopter_statefeedback_M->Timing.clockTickH0 *
        helicopter_statefeedback_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter_statefeedback_M)) {
    helicopter_statefeedback_M->Timing.t[0] = rtsiGetT
      (&helicopter_statefeedback_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter_statefeedback/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicopter_statefeedback_DW.HILReadEncoderTimebase_Task, 1,
         &helicopter_statefeedback_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter_statefeedback_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter_statefeedback_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter_statefeedback_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S4>/Travel: Count to rad' */
    helicopter_statefeedback_B.TravelCounttorad =
      helicopter_statefeedback_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S11>/Gain' */
    helicopter_statefeedback_B.Gain = helicopter_statefeedback_P.Gain_Gain *
      helicopter_statefeedback_B.TravelCounttorad;

    /* Sum: '<Root>/Sum1' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]1'
     */
    helicopter_statefeedback_B.Sum1 =
      helicopter_statefeedback_P.elavation_offsetdeg1_Value +
      helicopter_statefeedback_B.Gain;
  }

  /* TransferFcn: '<S4>/Travel: Transfer Fcn' */
  rtb_Derivative = 0.0;
  rtb_Derivative += helicopter_statefeedback_P.TravelTransferFcn_C *
    helicopter_statefeedback_X.TravelTransferFcn_CSTATE;
  rtb_Derivative += helicopter_statefeedback_P.TravelTransferFcn_D *
    helicopter_statefeedback_B.TravelCounttorad;

  /* Gain: '<S12>/Gain' */
  helicopter_statefeedback_B.Gain_d = helicopter_statefeedback_P.Gain_Gain_l *
    rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter_statefeedback_B.PitchCounttorad =
      helicopter_statefeedback_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S8>/Gain' */
    helicopter_statefeedback_B.Gain_i = helicopter_statefeedback_P.Gain_Gain_a *
      helicopter_statefeedback_B.PitchCounttorad;
  }

  /* TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  rtb_Derivative = 0.0;
  rtb_Derivative += helicopter_statefeedback_P.PitchTransferFcn_C *
    helicopter_statefeedback_X.PitchTransferFcn_CSTATE;
  rtb_Derivative += helicopter_statefeedback_P.PitchTransferFcn_D *
    helicopter_statefeedback_B.PitchCounttorad;

  /* Gain: '<S9>/Gain' */
  helicopter_statefeedback_B.Gain_b = helicopter_statefeedback_P.Gain_Gain_ae *
    rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    helicopter_statefeedback_B.ElevationCounttorad =
      helicopter_statefeedback_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S6>/Gain' */
    helicopter_statefeedback_B.Gain_e = helicopter_statefeedback_P.Gain_Gain_lv *
      helicopter_statefeedback_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter_statefeedback_B.Sum = helicopter_statefeedback_B.Gain_e +
      helicopter_statefeedback_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  rtb_Derivative = 0.0;
  rtb_Derivative += helicopter_statefeedback_P.ElevationTransferFcn_C *
    helicopter_statefeedback_X.ElevationTransferFcn_CSTATE;
  rtb_Derivative += helicopter_statefeedback_P.ElevationTransferFcn_D *
    helicopter_statefeedback_B.ElevationCounttorad;

  /* Gain: '<S7>/Gain' */
  helicopter_statefeedback_B.Gain_dg = helicopter_statefeedback_P.Gain_Gain_n *
    rtb_Derivative;

  /* Gain: '<S2>/Gain1' */
  helicopter_statefeedback_B.Gain1[0] = helicopter_statefeedback_P.Gain1_Gain *
    helicopter_statefeedback_B.Sum1;
  helicopter_statefeedback_B.Gain1[1] = helicopter_statefeedback_P.Gain1_Gain *
    helicopter_statefeedback_B.Gain_d;
  helicopter_statefeedback_B.Gain1[2] = helicopter_statefeedback_P.Gain1_Gain *
    helicopter_statefeedback_B.Gain_i;
  helicopter_statefeedback_B.Gain1[3] = helicopter_statefeedback_P.Gain1_Gain *
    helicopter_statefeedback_B.Gain_b;
  helicopter_statefeedback_B.Gain1[4] = helicopter_statefeedback_P.Gain1_Gain *
    helicopter_statefeedback_B.Sum;
  helicopter_statefeedback_B.Gain1[5] = helicopter_statefeedback_P.Gain1_Gain *
    helicopter_statefeedback_B.Gain_dg;
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter_statefeedback_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter_statefeedback_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex =
      helicopter_statefeedback_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter_statefeedback_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    helicopter_statefeedback_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Derivative = pDataValues[currTimeIndex];
        } else {
          rtb_Derivative = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Derivative = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helicopter_statefeedback_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter_statefeedback_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex =
      helicopter_statefeedback_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helicopter_statefeedback_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    helicopter_statefeedback_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&rtb_Sum2_f[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&rtb_Sum2_f[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 141;
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
          for (elIdx = 0; elIdx < 4; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum2_f[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  /* Sum: '<Root>/Sum2' */
  rtb_Sum2_f[0] = helicopter_statefeedback_B.Gain1[0] - rtb_Sum2_f[0];
  rtb_Sum2_f[1] = helicopter_statefeedback_B.Gain1[1] - rtb_Sum2_f[1];
  rtb_Sum2_f[2] = helicopter_statefeedback_B.Gain1[2] - rtb_Sum2_f[2];
  rtb_Sum2_f[3] = helicopter_statefeedback_B.Gain1[3] - rtb_Sum2_f[3];

  /* Gain: '<Root>/Gain' */
  rtb_Backgain = ((helicopter_statefeedback_P.K[0] * rtb_Sum2_f[0] +
                   helicopter_statefeedback_P.K[1] * rtb_Sum2_f[1]) +
                  helicopter_statefeedback_P.K[2] * rtb_Sum2_f[2]) +
    helicopter_statefeedback_P.K[3] * rtb_Sum2_f[3];

  /* Sum: '<S5>/Sum' incorporates:
   *  Constant: '<S5>/Vd_bias'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<Root>/Sum3'
   *  Sum: '<S5>/Sum2'
   *  Sum: '<S5>/Sum3'
   */
  helicopter_statefeedback_B.Sum_k = (((rtb_Derivative - rtb_Backgain) -
    helicopter_statefeedback_B.Gain1[2]) * helicopter_statefeedback_P.K_pp -
    helicopter_statefeedback_P.K_pd * helicopter_statefeedback_B.Gain1[3]) +
    helicopter_statefeedback_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
  }

  /* Integrator: '<S3>/Integrator'
   *
   * Regarding '<S3>/Integrator':
   *  Limited Integrator
   */
  if (helicopter_statefeedback_X.Integrator_CSTATE >=
      helicopter_statefeedback_P.Integrator_UpperSat ) {
    helicopter_statefeedback_X.Integrator_CSTATE =
      helicopter_statefeedback_P.Integrator_UpperSat;
  } else if (helicopter_statefeedback_X.Integrator_CSTATE <=
             (helicopter_statefeedback_P.Integrator_LowerSat) ) {
    helicopter_statefeedback_X.Integrator_CSTATE =
      (helicopter_statefeedback_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter_statefeedback_X.Integrator_CSTATE;

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   */
  rtb_Derivative = helicopter_statefeedback_P.elevation_ref_Value -
    helicopter_statefeedback_B.Gain1[4];

  /* Sum: '<S3>/Sum2' incorporates:
   *  Constant: '<S3>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  helicopter_statefeedback_B.Sum2 = ((helicopter_statefeedback_P.K_ep *
    rtb_Derivative + rtb_Backgain) - helicopter_statefeedback_P.K_ed *
    helicopter_statefeedback_B.Gain1[5]) + helicopter_statefeedback_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helicopter_statefeedback_B.Sum2 -
                  helicopter_statefeedback_B.Sum_k) *
    helicopter_statefeedback_P.Backgain_Gain;

  /* Gain: '<S3>/K_ei' */
  helicopter_statefeedback_B.K_ei = helicopter_statefeedback_P.K_ei *
    rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helicopter_statefeedback_DW.TimeStampA >=
       helicopter_statefeedback_M->Timing.t[0]) &&
      (helicopter_statefeedback_DW.TimeStampB >=
       helicopter_statefeedback_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    lastTime = helicopter_statefeedback_DW.TimeStampA;
    lastU = &helicopter_statefeedback_DW.LastUAtTimeA;
    if (helicopter_statefeedback_DW.TimeStampA <
        helicopter_statefeedback_DW.TimeStampB) {
      if (helicopter_statefeedback_DW.TimeStampB <
          helicopter_statefeedback_M->Timing.t[0]) {
        lastTime = helicopter_statefeedback_DW.TimeStampB;
        lastU = &helicopter_statefeedback_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter_statefeedback_DW.TimeStampA >=
          helicopter_statefeedback_M->Timing.t[0]) {
        lastTime = helicopter_statefeedback_DW.TimeStampB;
        lastU = &helicopter_statefeedback_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter_statefeedback_B.PitchCounttorad - *lastU) /
      (helicopter_statefeedback_M->Timing.t[0] - lastTime);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S10>/Gain' */
  helicopter_statefeedback_B.Gain_l = helicopter_statefeedback_P.Gain_Gain_a1 *
    rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
  }

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter_statefeedback_P.BackmotorSaturation_UpperSat) {
    helicopter_statefeedback_B.BackmotorSaturation =
      helicopter_statefeedback_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain <
             helicopter_statefeedback_P.BackmotorSaturation_LowerSat) {
    helicopter_statefeedback_B.BackmotorSaturation =
      helicopter_statefeedback_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter_statefeedback_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  lastTime = (helicopter_statefeedback_B.Sum_k + helicopter_statefeedback_B.Sum2)
    * helicopter_statefeedback_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (lastTime > helicopter_statefeedback_P.FrontmotorSaturation_UpperSat) {
    helicopter_statefeedback_B.FrontmotorSaturation =
      helicopter_statefeedback_P.FrontmotorSaturation_UpperSat;
  } else if (lastTime < helicopter_statefeedback_P.FrontmotorSaturation_LowerSat)
  {
    helicopter_statefeedback_B.FrontmotorSaturation =
      helicopter_statefeedback_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter_statefeedback_B.FrontmotorSaturation = lastTime;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter_statefeedback/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter_statefeedback_DW.HILWriteAnalog_Buffer[0] =
        helicopter_statefeedback_B.FrontmotorSaturation;
      helicopter_statefeedback_DW.HILWriteAnalog_Buffer[1] =
        helicopter_statefeedback_B.BackmotorSaturation;
      result = hil_write_analog(helicopter_statefeedback_DW.HILInitialize_Card,
        helicopter_statefeedback_P.HILWriteAnalog_channels, 2,
        &helicopter_statefeedback_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter_statefeedback_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter_statefeedback_DW.TimeStampA == (rtInf)) {
    helicopter_statefeedback_DW.TimeStampA =
      helicopter_statefeedback_M->Timing.t[0];
    lastU = &helicopter_statefeedback_DW.LastUAtTimeA;
  } else if (helicopter_statefeedback_DW.TimeStampB == (rtInf)) {
    helicopter_statefeedback_DW.TimeStampB =
      helicopter_statefeedback_M->Timing.t[0];
    lastU = &helicopter_statefeedback_DW.LastUAtTimeB;
  } else if (helicopter_statefeedback_DW.TimeStampA <
             helicopter_statefeedback_DW.TimeStampB) {
    helicopter_statefeedback_DW.TimeStampA =
      helicopter_statefeedback_M->Timing.t[0];
    lastU = &helicopter_statefeedback_DW.LastUAtTimeA;
  } else {
    helicopter_statefeedback_DW.TimeStampB =
      helicopter_statefeedback_M->Timing.t[0];
    lastU = &helicopter_statefeedback_DW.LastUAtTimeB;
  }

  *lastU = helicopter_statefeedback_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter_statefeedback_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter_statefeedback_M->solverInfo);
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
  if (!(++helicopter_statefeedback_M->Timing.clockTick0)) {
    ++helicopter_statefeedback_M->Timing.clockTickH0;
  }

  helicopter_statefeedback_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter_statefeedback_M->solverInfo);

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
    if (!(++helicopter_statefeedback_M->Timing.clockTick1)) {
      ++helicopter_statefeedback_M->Timing.clockTickH1;
    }

    helicopter_statefeedback_M->Timing.t[1] =
      helicopter_statefeedback_M->Timing.clockTick1 *
      helicopter_statefeedback_M->Timing.stepSize1 +
      helicopter_statefeedback_M->Timing.clockTickH1 *
      helicopter_statefeedback_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter_statefeedback_derivatives(void)
{
  XDot_helicopter_statefeedback_T *_rtXdot;
  _rtXdot = ((XDot_helicopter_statefeedback_T *)
             helicopter_statefeedback_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE +=
    helicopter_statefeedback_P.TravelTransferFcn_A *
    helicopter_statefeedback_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE +=
    helicopter_statefeedback_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE +=
    helicopter_statefeedback_P.PitchTransferFcn_A *
    helicopter_statefeedback_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_statefeedback_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter_statefeedback_P.ElevationTransferFcn_A *
    helicopter_statefeedback_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter_statefeedback_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter_statefeedback_X.Integrator_CSTATE <=
            (helicopter_statefeedback_P.Integrator_LowerSat) );
    usat = ( helicopter_statefeedback_X.Integrator_CSTATE >=
            helicopter_statefeedback_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter_statefeedback_B.K_ei > 0)) ||
        (usat && (helicopter_statefeedback_B.K_ei < 0)) ) {
      ((XDot_helicopter_statefeedback_T *)
        helicopter_statefeedback_M->ModelData.derivs)->Integrator_CSTATE =
        helicopter_statefeedback_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter_statefeedback_T *)
        helicopter_statefeedback_M->ModelData.derivs)->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter_statefeedback_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_statefeedback/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0",
                      &helicopter_statefeedback_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (helicopter_statefeedback_DW.HILInitialize_Card,
       "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter_statefeedback_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
      return;
    }

    if ((helicopter_statefeedback_P.HILInitialize_set_analog_input_ &&
         !is_switching) ||
        (helicopter_statefeedback_P.HILInitialize_set_analog_inpu_m &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums =
          &helicopter_statefeedback_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] =
            helicopter_statefeedback_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums =
          &helicopter_statefeedback_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] =
            helicopter_statefeedback_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges
        (helicopter_statefeedback_DW.HILInitialize_Card,
         helicopter_statefeedback_P.HILInitialize_analog_input_chan, 8U,
         &helicopter_statefeedback_DW.HILInitialize_AIMinimums[0],
         &helicopter_statefeedback_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_statefeedback_P.HILInitialize_set_analog_output &&
         !is_switching) ||
        (helicopter_statefeedback_P.HILInitialize_set_analog_outp_b &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums =
          &helicopter_statefeedback_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] =
            helicopter_statefeedback_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums =
          &helicopter_statefeedback_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] =
            helicopter_statefeedback_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges
        (helicopter_statefeedback_DW.HILInitialize_Card,
         helicopter_statefeedback_P.HILInitialize_analog_output_cha, 8U,
         &helicopter_statefeedback_DW.HILInitialize_AOMinimums[0],
         &helicopter_statefeedback_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_statefeedback_P.HILInitialize_set_analog_outp_e &&
         !is_switching) ||
        (helicopter_statefeedback_P.HILInitialize_set_analog_outp_j &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &helicopter_statefeedback_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter_statefeedback_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter_statefeedback_DW.HILInitialize_Card,
        helicopter_statefeedback_P.HILInitialize_analog_output_cha, 8U,
        &helicopter_statefeedback_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_statefeedback_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &helicopter_statefeedback_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter_statefeedback_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter_statefeedback_DW.HILInitialize_Card,
         helicopter_statefeedback_P.HILInitialize_analog_output_cha, 8U,
         &helicopter_statefeedback_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_statefeedback_P.HILInitialize_set_encoder_param &&
         !is_switching) ||
        (helicopter_statefeedback_P.HILInitialize_set_encoder_par_m &&
         is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter_statefeedback_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] =
            helicopter_statefeedback_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter_statefeedback_DW.HILInitialize_Card,
         helicopter_statefeedback_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter_statefeedback_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_statefeedback_P.HILInitialize_set_encoder_count &&
         !is_switching) ||
        (helicopter_statefeedback_P.HILInitialize_set_encoder_cou_k &&
         is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter_statefeedback_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter_statefeedback_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts
        (helicopter_statefeedback_DW.HILInitialize_Card,
         helicopter_statefeedback_P.HILInitialize_encoder_channels, 8U,
         &helicopter_statefeedback_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_statefeedback_P.HILInitialize_set_pwm_params_at &&
         !is_switching) ||
        (helicopter_statefeedback_P.HILInitialize_set_pwm_params__f &&
         is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_statefeedback_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            helicopter_statefeedback_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter_statefeedback_DW.HILInitialize_Card,
        helicopter_statefeedback_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter_statefeedback_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter_statefeedback_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues =
          &helicopter_statefeedback_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter_statefeedback_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter_statefeedback_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = helicopter_statefeedback_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter_statefeedback_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = p_HILInitialize_pwm_channels[i1];
            helicopter_statefeedback_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              helicopter_statefeedback_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency
          (helicopter_statefeedback_DW.HILInitialize_Card,
           &helicopter_statefeedback_DW.HILInitialize_POSortedChans[0],
           num_duty_cycle_modes,
           &helicopter_statefeedback_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle
          (helicopter_statefeedback_DW.HILInitialize_Card,
           &helicopter_statefeedback_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
           num_frequency_modes,
           &helicopter_statefeedback_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_statefeedback_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            helicopter_statefeedback_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter_statefeedback_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] =
            helicopter_statefeedback_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter_statefeedback_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] =
            helicopter_statefeedback_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration
        (helicopter_statefeedback_DW.HILInitialize_Card,
         helicopter_statefeedback_P.HILInitialize_pwm_channels, 8U,
         (t_pwm_configuration *)
         &helicopter_statefeedback_DW.HILInitialize_POModeValues[0],
         (t_pwm_alignment *)
         &helicopter_statefeedback_DW.HILInitialize_POAlignValues[0],
         (t_pwm_polarity *)
         &helicopter_statefeedback_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicopter_statefeedback_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] =
            helicopter_statefeedback_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues =
          &helicopter_statefeedback_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            helicopter_statefeedback_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband
        (helicopter_statefeedback_DW.HILInitialize_Card,
         helicopter_statefeedback_P.HILInitialize_pwm_channels, 8U,
         &helicopter_statefeedback_DW.HILInitialize_POSortedFreqs[0],
         &helicopter_statefeedback_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_statefeedback_P.HILInitialize_set_pwm_outputs_a &&
         !is_switching) ||
        (helicopter_statefeedback_P.HILInitialize_set_pwm_outputs_g &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues =
          &helicopter_statefeedback_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            helicopter_statefeedback_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter_statefeedback_DW.HILInitialize_Card,
        helicopter_statefeedback_P.HILInitialize_pwm_channels, 8U,
        &helicopter_statefeedback_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_statefeedback_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues =
          &helicopter_statefeedback_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            helicopter_statefeedback_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter_statefeedback_DW.HILInitialize_Card,
         helicopter_statefeedback_P.HILInitialize_pwm_channels, 8U,
         &helicopter_statefeedback_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter_statefeedback/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (helicopter_statefeedback_DW.HILInitialize_Card,
       helicopter_statefeedback_P.HILReadEncoderTimebase_samples_,
       helicopter_statefeedback_P.HILReadEncoderTimebase_channels, 3,
       &helicopter_statefeedback_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
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
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559807444,
      0.52359877559803447, 0.52359877559798307, 0.52359877559791523,
      0.523598775597823, 0.5235987755976923, 0.52359877559749723,
      0.52359877559718215, 0.5235987755966055, 0.52359877559526735,
      0.52359877558925361, 0.38860546189573791, 0.10951698887526917,
      -0.1100383394610223, -0.27691123755291591, -0.3979068273164692,
      -0.47963675127901134, -0.523591909441963, -0.52359877534342569,
      -0.52359877536343258, -0.5235987224145312, -0.50343489391626584,
      -0.46497044043306424, -0.42077345861352544, -0.37347862555446432,
      -0.32522986550994049, -0.27772344682220712, -0.2322553027667397,
      -0.1897701881732447, -0.15091074088283754, -0.11606494657763659,
      -0.08541089068498342, -0.058958016011036223, -0.036584388943327412,
      -0.018069712651425172, -0.0031240162861197467, 0.008587900985857537,
      0.017426078578606438, 0.023758780967378875, 0.027949332008985642,
      0.030346034333594667, 0.03127479682777045, 0.0310340938210644,
      0.0298918916543043, 0.02808419891071846, 0.025814923189764548,
      0.023256747740391905, 0.020552773751695553, 0.017818707168261894,
      0.015145401405518497, 0.01260159841128787, 0.010236739517233149,
      0.0080837440173715247, 0.0061616771427812711, 0.0044782499569109863,
      0.0030321116735972076, 0.0018149100878760917, 0.00081310836210016267,
      9.5565276842160941E-6, -0.00061517602288165755, -0.0010816948208518151,
      -0.0014108848952404636, -0.0016232063484833872, -0.0017381600630382574,
      -0.0017739020793310287, -0.0017469853873323285, -0.0016722086841820883,
      -0.0015625529081204723, -0.0014291879259837978, -0.0012815335105175723,
      -0.0011273605978063502, -0.00097292068600022858, -0.00082309306297666392,
      -0.00068154128644839379, -0.00055087195231493555, -0.00043279025403042593,
      -0.00032824814515925653, -0.00023758206462032373, -0.00016063817112712349,
      -9.688386664573433E-5, -4.5505078677972685E-5, -5.4893309923764594E-6,
      2.4304922656237162E-5, 4.5091880270843108E-5, 5.8113799538615989E-5,
      6.4606933573472039E-5, 6.5776522273753571E-5, 6.2779201128075149E-5,
      5.6710998483779331E-5, 4.8598837981454231E-5, 3.9393103171466183E-5,
      2.9958323012006171E-5, 2.1058409934228412E-5, 1.3332266149043419E-5,
      7.2554391803702282E-6, 3.0850995699961848E-6, 7.9190498794152739E-7,
      -4.0593744014257015E-17, -1.5332144970613159E-17, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

    helicopter_statefeedback_DW.FromWorkspace_PWORK.TimePtr = (void *)
      pTimeValues0;
    helicopter_statefeedback_DW.FromWorkspace_PWORK.DataPtr = (void *)
      pDataValues0;
    helicopter_statefeedback_DW.FromWorkspace_IWORK.PrevIndex = 0;
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
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413625279, 3.1262155534580049,
      3.1033093000299821, 3.0666274151912156, 3.014453922394225,
      2.9456562771176764, 2.8595077632937147, 2.7555515879654058,
      2.6335051104906526, 2.4931956060326308, 2.3345185760651095,
      2.1583782719260478, 1.9678000778704727, 1.7674110445466207,
      1.5625810512761342, 1.3586839679386826, 1.1605920903248919,
      0.97235419679771273, 0.79695012325984049, 0.63643293374959042,
      0.492159597969274, 0.36485729820996032, 0.25463993267206,
      0.1610962359301791, 0.083400591252193013, 0.020423624533390517,
      -0.029167512815941744, -0.066823004653064969, -0.094039185304793951,
      -0.1123016902197959, -0.12304167596519094, -0.12760358353335188,
      -0.12722285379205739, -0.12301203948381061, -0.11595384747675903,
      -0.10689976340781571, -0.096573044536239716, -0.0855750068726116,
      -0.074393673542786518, -0.0634139886026101, -0.052928930978585753,
      -0.0431509845929514, -0.034223531491125975, -0.026231833998759255,
      -0.019213359191285469, -0.013167274288011305, -0.0080630053551753473,
      -0.0038478045493674498, -0.00045331387359379306, 0.0021988530036483648,
      0.0041924638884134058, 0.0056128908698021675, 0.0065441252851624956,
      0.0070665002189294682, 0.007255032990408989, 0.00717830051005529,
      0.006897763395254222, 0.0064674596661119165, 0.00593399511017771,
      0.0053367645211473752, 0.0047083455725129373, 0.0040750147523033915,
      0.0034573422992259209, 0.0028708302460222226, 0.0023265643498006218,
      0.0018318567755375491, 0.0013908618413420552, 0.0010051519083085054,
      0.00067424460561774189, 0.00039607604483589968, 0.00016741753288223544,
      -1.5764411858470837E-5, -0.00015800314738642176, -0.00026407916771299997,
      -0.00033881565345682897, -0.00038691781724583353, -0.0004128512108077269,
      -0.00042075405089002696, -0.00041437869183475569, -0.00039705758179327467,
      -0.00037168935026994572, -0.0003407410549825221, -0.0003062630387661212,
      -0.00026991328973278837, -0.00023298864146416003, -0.00019646057949526021,
      -0.00016101382346239446, -0.00012708622092749209, -9.4908810282185864E-5,
      -6.4545177863868E-5, -3.5929439305405029E-5, -8.90230607054186E-6,
      1.6755258848342884E-5, 4.12913533646134E-5, 6.4956154500282433E-5,
      8.7979679297890866E-5, 0.00011055569399486016, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048909061991,
      -0.046506351618091031, -0.091625013712091086, -0.14672753935506733,
      -0.20869397118796093, -0.2751905811061946, -0.34459405529584686,
      -0.415824701313235, -0.48818590989901212, -0.5612380178320886,
      -0.63470811987008435, -0.70456121655624571, -0.7623127762223002,
      -0.80155613329540887, -0.81931997308194582, -0.81558833334980685,
      -0.79236751045516307, -0.75295157410871649, -0.70161629415148918,
      -0.642068758041, -0.5770933431212657, -0.50920919903725481,
      -0.4408694621516015, -0.37417478696752354, -0.31078257871194431,
      -0.25190786687520994, -0.19836454939732906, -0.1506219673484929,
      -0.10886472260691592, -0.073050019660007848, -0.042959942981580151,
      -0.018247630272643797, 0.0015229189651780261, 0.016843257232987102,
      0.028232768028206262, 0.036216336275773323, 0.041306875486303914,
      0.043992150654512473, 0.044725333319300294, 0.043918739760705684,
      0.041940230496097373, 0.039111785542537428, 0.0357098124073017,
      0.031966789969466884, 0.028073899229895151, 0.024184339613096656,
      0.020417075731343829, 0.016860803223231587, 0.013577962703094627,
      0.010608667508968632, 0.007974443539060164, 0.0056817079255550469,
      0.0037249376614413112, 0.0020894997350678905, 0.00075413108591808223,
      -0.00030692992141479565, -0.00112214845920427, -0.0017212149165692226,
      -0.002133858223736824, -0.0023889223561213407, -0.0025136757945377537,
      -0.0025333232808381809, -0.0024706898123098837, -0.0023460482128147947,
      -0.0021770635848864044, -0.0019788302970522903, -0.0017639797367819756,
      -0.0015428397321341988, -0.0013236292107630541, -0.0011126742431273691,
      -0.000914634047814657, -0.00073272777896282509, -0.00056895494211180371,
      -0.00042430408130631283, -0.00029894594297531606, -0.00019240865515601812,
      -0.00010373357424757337, -3.1611360329200274E-5, 2.550143622108505E-5,
      6.9284440165924E-5, 0.0001014729260933157, 0.00012379318114969477,
      0.00013791206486560344, 0.00014539899613333143, 0.00014769859307451336,
      0.00014611224787559922, 0.000141787024131463, 0.00013571041013960949,
      0.00012870964258122485, 0.00012145452967327145, 0.0001144629542338519,
      0.00010810853293945266, 0.00010263025967553898, 9.8144378065082063E-5,
      9.465920454267608E-5, 9.2094099190433756E-5, 9.0304058787877176E-5,
      8.91109587141194E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.10602875205861005, 0.22266037932307303, 0.31888147181624232,
      0.38944360631121483, 0.43795507377648213, 0.46997264230352054,
      0.49051724877497993, 0.50343100141409236, 0.51142138585938279,
      0.51630439857559973, 0.5192586212675202, 0.49369500885904027,
      0.4081659670587614, 0.27735705984392089, 0.12554803518857219,
      -0.02637380442694166, -0.16411590764817871, -0.2785767842358442,
      -0.36281815260289824, -0.42085924264294028, -0.45922141703380015,
      -0.47977920385573453, -0.482999159776279, -0.47137249195678604,
      -0.44803191699491934, -0.41610397764324641, -0.37842371849810491,
      -0.33742633592109822, -0.29512425777261336, -0.25312464196251705,
      -0.21266516986452466, -0.17465718802102487, -0.13973069118134349,
      -0.10827829968414161, -0.080496712382112978, -0.05642481126381068,
      -0.035977986944966428, -0.018978499319992535, -0.0051818550552452224,
      0.0057006952153475883, 0.013983347843646946, 0.019990368682992511,
      0.024043846827354103, 0.026454253043480896, 0.027513464961974859,
      0.027489921858446696, 0.026625582310723334, 0.02513437586433389,
      0.023201863003346503, 0.020985844389337649, 0.018617688948044889,
      0.016204179667314456, 0.013829704890780056, 0.011558650651908645,
      0.0094378756039130169, 0.007499173956005389, 0.0057616532742171484,
      0.0042339729232725165, 0.002916405297005134, 0.0018026958727821757,
      0.00088170965649384762, 0.00013886092934043047, -0.00044266943438280129,
      -0.00088091922171180872, -0.0011943188110467254, -0.0014010371685183724,
      -0.001518481703575448, -0.0015629330943505808, -0.0015492962434657848,
      -0.001490949142651121, -0.0013996724643236191, -0.0012856440340190632,
      -0.001157483862215849, -0.001022337038652467, -0.00088598344454119391,
      -0.00075296486124410322, -0.00062672160478955178, -0.00050973226282778129,
      -0.00040365143331496246, -0.00030944155013681613, -0.000227495925005011,
      -0.00015775103811508149, -9.9786877774685792E-5, -5.2914770767496487E-5,
      -1.625267291083395E-5, 1.1211682003878427E-5, 3.0569029520228633E-5,
      4.294718689451349E-5, 4.9478751347060993E-5, 5.1276367137421242E-5,
      4.9413784960899323E-5, 4.4910622807850582E-5, 3.8718343149644846E-5,
      3.1704497960415828E-5, 2.4631875387714816E-5, 1.8129184956384737E-5,
      1.265132190729764E-5, 8.4323756487128E-6, 5.4485883366234407E-6, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823444021,
      0.46652650905785192, 0.38488436997267728, 0.28224853797988991,
      0.19404586986106923, 0.12807027410815364, 0.082178425885837569,
      0.051655010556449679, 0.031961537781161606, 0.019532050864867829,
      0.011816890767682156, -0.10225444963391987, -0.34211616720111526,
      -0.52323562885936215, -0.60723609862139483, -0.60768735846205535,
      -0.55096841288494813, -0.457843506350662, -0.33696547346821631,
      -0.23216436016016812, -0.15344869756343951, -0.0822311472877376,
      -0.012879823682177707, 0.046506671277971548, 0.093362299847467,
      0.12771175740669166, 0.15072103658056601, 0.16398953030802704,
      0.16920831259393943, 0.16799846324038506, 0.16183788839196964,
      0.1520319273739992, 0.13970598735872553, 0.12580956598880746,
      0.11112634920811458, 0.096287604473209165, 0.081787297275377022,
      0.067997950499895557, 0.055186577058989245, 0.043530201082371239,
      0.033130610513197434, 0.024028083357382258, 0.016213912577446362,
      0.00964162486450718, 0.0042368476739758475, -9.4172414112656268E-5,
      -0.0034573581908934571, -0.00596482578555777, -0.0077300514439495429,
      -0.0088640744560354225, -0.009472621765171033, -0.0096540371229217464,
      -0.0094978991061375976, -0.0090842169554856377, -0.00848310019198251,
      -0.0077548065916305158, -0.0069500827271529609, -0.0061107214037785284,
      -0.0052702705050695284, -0.004454837696891833, -0.003683944865153312,
      -0.0029713949086136687, -0.0023261214548929267, -0.0017529991493160297,
      -0.0012535983573396671, -0.00082687342988658769, -0.00046977814022830278,
      -0.00017780556310053122, 5.4547403539183512E-5, 0.00023338840325865535,
      0.00036510671331000744, 0.00045611372121822361, 0.00051264068721285582,
      0.00054058729425352794, 0.00054541437644509269, 0.00053207433318836254,
      0.00050497302581820557, 0.00046795736784708186, 0.00042432331805127509,
      0.00037683953271258548, 0.0003277825005272203, 0.00027897954755971816,
      0.00023185664136158283, 0.00018748842802875722, 0.00014664839142665014,
      0.0001098574196588495, 7.742939006540083E-5, 4.9512629497139429E-5,
      2.6126257810190009E-5, 7.190463161441E-6, -7.4503287060876995E-6,
      -1.8012648612194976E-5, -2.4769118632822928E-5, -2.8055380756916086E-5,
      -2.8290490290804048E-5, -2.6010761725320315E-5, -2.1911452196348394E-5,
      -1.6875785034339357E-5, -1.1935149248357441E-5, -8.023739200293113E-6, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0 } ;

    helicopter_statefeedback_DW.FromWorkspace1_PWORK.TimePtr = (void *)
      pTimeValues0;
    helicopter_statefeedback_DW.FromWorkspace1_PWORK.DataPtr = (void *)
      pDataValues0;
    helicopter_statefeedback_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter_statefeedback_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter_statefeedback_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter_statefeedback_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter_statefeedback_X.Integrator_CSTATE =
    helicopter_statefeedback_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter_statefeedback_DW.TimeStampA = (rtInf);
  helicopter_statefeedback_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter_statefeedback_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_statefeedback/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter_statefeedback_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter_statefeedback_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter_statefeedback_P.HILInitialize_set_analog_out_ex &&
         !is_switching) ||
        (helicopter_statefeedback_P.HILInitialize_set_analog_outp_c &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &helicopter_statefeedback_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter_statefeedback_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter_statefeedback_P.HILInitialize_set_pwm_output_ap &&
         !is_switching) ||
        (helicopter_statefeedback_P.HILInitialize_set_pwm_outputs_p &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues =
          &helicopter_statefeedback_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            helicopter_statefeedback_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter_statefeedback_DW.HILInitialize_Card
                         ,
                         helicopter_statefeedback_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter_statefeedback_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         ,
                         &helicopter_statefeedback_DW.HILInitialize_AOVoltages[0]
                         , &helicopter_statefeedback_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (helicopter_statefeedback_DW.HILInitialize_Card,
             helicopter_statefeedback_P.HILInitialize_analog_output_cha,
             num_final_analog_outputs,
             &helicopter_statefeedback_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm
            (helicopter_statefeedback_DW.HILInitialize_Card,
             helicopter_statefeedback_P.HILInitialize_pwm_channels,
             num_final_pwm_outputs,
             &helicopter_statefeedback_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_statefeedback_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter_statefeedback_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter_statefeedback_DW.HILInitialize_Card);
    hil_close(helicopter_statefeedback_DW.HILInitialize_Card);
    helicopter_statefeedback_DW.HILInitialize_Card = NULL;
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
  helicopter_statefeedback_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter_statefeedback_update();
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
  helicopter_statefeedback_initialize();
}

void MdlTerminate(void)
{
  helicopter_statefeedback_terminate();
}

/* Registration function */
RT_MODEL_helicopter_statefeed_T *helicopter_statefeedback(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter_statefeedback_P.Integrator_UpperSat = rtInf;
  helicopter_statefeedback_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter_statefeedback_M, 0,
                sizeof(RT_MODEL_helicopter_statefeed_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter_statefeedback_M->solverInfo,
                          &helicopter_statefeedback_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter_statefeedback_M->solverInfo, &rtmGetTPtr
                (helicopter_statefeedback_M));
    rtsiSetStepSizePtr(&helicopter_statefeedback_M->solverInfo,
                       &helicopter_statefeedback_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter_statefeedback_M->solverInfo,
                 &helicopter_statefeedback_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter_statefeedback_M->solverInfo, (real_T **)
                         &helicopter_statefeedback_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter_statefeedback_M->solverInfo,
      &helicopter_statefeedback_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter_statefeedback_M->solverInfo,
                          (&rtmGetErrorStatus(helicopter_statefeedback_M)));
    rtsiSetRTModelPtr(&helicopter_statefeedback_M->solverInfo,
                      helicopter_statefeedback_M);
  }

  rtsiSetSimTimeStep(&helicopter_statefeedback_M->solverInfo, MAJOR_TIME_STEP);
  helicopter_statefeedback_M->ModelData.intgData.f[0] =
    helicopter_statefeedback_M->ModelData.odeF[0];
  helicopter_statefeedback_M->ModelData.contStates = ((real_T *)
    &helicopter_statefeedback_X);
  rtsiSetSolverData(&helicopter_statefeedback_M->solverInfo, (void *)
                    &helicopter_statefeedback_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter_statefeedback_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter_statefeedback_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter_statefeedback_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter_statefeedback_M->Timing.sampleTimes =
      (&helicopter_statefeedback_M->Timing.sampleTimesArray[0]);
    helicopter_statefeedback_M->Timing.offsetTimes =
      (&helicopter_statefeedback_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter_statefeedback_M->Timing.sampleTimes[0] = (0.0);
    helicopter_statefeedback_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter_statefeedback_M->Timing.offsetTimes[0] = (0.0);
    helicopter_statefeedback_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter_statefeedback_M,
             &helicopter_statefeedback_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter_statefeedback_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter_statefeedback_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter_statefeedback_M, 35.0);
  helicopter_statefeedback_M->Timing.stepSize0 = 0.002;
  helicopter_statefeedback_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter_statefeedback_M->Sizes.checksums[0] = (3607639473U);
  helicopter_statefeedback_M->Sizes.checksums[1] = (4003152393U);
  helicopter_statefeedback_M->Sizes.checksums[2] = (482794593U);
  helicopter_statefeedback_M->Sizes.checksums[3] = (546085011U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter_statefeedback_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter_statefeedback_M->extModeInfo,
      &helicopter_statefeedback_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter_statefeedback_M->extModeInfo,
                        helicopter_statefeedback_M->Sizes.checksums);
    rteiSetTPtr(helicopter_statefeedback_M->extModeInfo, rtmGetTPtr
                (helicopter_statefeedback_M));
  }

  helicopter_statefeedback_M->solverInfoPtr =
    (&helicopter_statefeedback_M->solverInfo);
  helicopter_statefeedback_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter_statefeedback_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter_statefeedback_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter_statefeedback_M->ModelData.blockIO = ((void *)
    &helicopter_statefeedback_B);

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      helicopter_statefeedback_B.Gain1[i] = 0.0;
    }

    helicopter_statefeedback_B.TravelCounttorad = 0.0;
    helicopter_statefeedback_B.Gain = 0.0;
    helicopter_statefeedback_B.Sum1 = 0.0;
    helicopter_statefeedback_B.Gain_d = 0.0;
    helicopter_statefeedback_B.PitchCounttorad = 0.0;
    helicopter_statefeedback_B.Gain_i = 0.0;
    helicopter_statefeedback_B.Gain_b = 0.0;
    helicopter_statefeedback_B.ElevationCounttorad = 0.0;
    helicopter_statefeedback_B.Gain_e = 0.0;
    helicopter_statefeedback_B.Sum = 0.0;
    helicopter_statefeedback_B.Gain_dg = 0.0;
    helicopter_statefeedback_B.Sum_k = 0.0;
    helicopter_statefeedback_B.Sum2 = 0.0;
    helicopter_statefeedback_B.K_ei = 0.0;
    helicopter_statefeedback_B.Gain_l = 0.0;
    helicopter_statefeedback_B.BackmotorSaturation = 0.0;
    helicopter_statefeedback_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter_statefeedback_M->ModelData.defaultParam = ((real_T *)
    &helicopter_statefeedback_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter_statefeedback_X;
    helicopter_statefeedback_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter_statefeedback_X, 0,
                  sizeof(X_helicopter_statefeedback_T));
  }

  /* states (dwork) */
  helicopter_statefeedback_M->ModelData.dwork = ((void *)
    &helicopter_statefeedback_DW);
  (void) memset((void *)&helicopter_statefeedback_DW, 0,
                sizeof(DW_helicopter_statefeedback_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_statefeedback_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_statefeedback_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_statefeedback_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_statefeedback_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_statefeedback_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_statefeedback_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_statefeedback_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_statefeedback_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter_statefeedback_DW.TimeStampA = 0.0;
  helicopter_statefeedback_DW.LastUAtTimeA = 0.0;
  helicopter_statefeedback_DW.TimeStampB = 0.0;
  helicopter_statefeedback_DW.LastUAtTimeB = 0.0;
  helicopter_statefeedback_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter_statefeedback_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter_statefeedback_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter_statefeedback_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter_statefeedback_M->Sizes.numY = (0);/* Number of model outputs */
  helicopter_statefeedback_M->Sizes.numU = (0);/* Number of model inputs */
  helicopter_statefeedback_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter_statefeedback_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter_statefeedback_M->Sizes.numBlocks = (60);/* Number of blocks */
  helicopter_statefeedback_M->Sizes.numBlockIO = (18);/* Number of block outputs */
  helicopter_statefeedback_M->Sizes.numBlockPrms = (146);/* Sum of parameter "widths" */
  return helicopter_statefeedback_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
