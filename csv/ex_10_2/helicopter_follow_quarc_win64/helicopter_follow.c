/*
 * helicopter_follow.c
 *
 * Code generation for model "helicopter_follow".
 *
 * Model version              : 1.168
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Fri Apr 21 14:38:15 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter_follow.h"
#include "helicopter_follow_private.h"
#include "helicopter_follow_dt.h"

/* Block signals (auto storage) */
B_helicopter_follow_T helicopter_follow_B;

/* Continuous states */
X_helicopter_follow_T helicopter_follow_X;

/* Block states (auto storage) */
DW_helicopter_follow_T helicopter_follow_DW;

/* Real-time model */
RT_MODEL_helicopter_follow_T helicopter_follow_M_;
RT_MODEL_helicopter_follow_T *const helicopter_follow_M = &helicopter_follow_M_;

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
  helicopter_follow_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter_follow_output(void)
{
  /* local block i/o variables */
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
    /* set solver stop time */
    if (!(helicopter_follow_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter_follow_M->solverInfo,
                            ((helicopter_follow_M->Timing.clockTickH0 + 1) *
        helicopter_follow_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter_follow_M->solverInfo,
                            ((helicopter_follow_M->Timing.clockTick0 + 1) *
        helicopter_follow_M->Timing.stepSize0 +
        helicopter_follow_M->Timing.clockTickH0 *
        helicopter_follow_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter_follow_M)) {
    helicopter_follow_M->Timing.t[0] = rtsiGetT(&helicopter_follow_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter_follow/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicopter_follow_DW.HILReadEncoderTimebase_Task, 1,
         &helicopter_follow_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter_follow_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter_follow_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter_follow_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S4>/Travel: Count to rad' */
    helicopter_follow_B.TravelCounttorad =
      helicopter_follow_P.TravelCounttorad_Gain * rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S11>/Gain' */
    helicopter_follow_B.Gain = helicopter_follow_P.Gain_Gain *
      helicopter_follow_B.TravelCounttorad;

    /* Sum: '<Root>/Sum1' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]1'
     */
    helicopter_follow_B.Sum1 = helicopter_follow_P.elavation_offsetdeg1_Value +
      helicopter_follow_B.Gain;
  }

  /* TransferFcn: '<S4>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter_follow_P.TravelTransferFcn_C *
    helicopter_follow_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helicopter_follow_P.TravelTransferFcn_D *
    helicopter_follow_B.TravelCounttorad;

  /* Gain: '<S12>/Gain' */
  helicopter_follow_B.Gain_d = helicopter_follow_P.Gain_Gain_l * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter_follow_B.PitchCounttorad =
      helicopter_follow_P.PitchCounttorad_Gain * rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S8>/Gain' */
    helicopter_follow_B.Gain_i = helicopter_follow_P.Gain_Gain_a *
      helicopter_follow_B.PitchCounttorad;
  }

  /* TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter_follow_P.PitchTransferFcn_C *
    helicopter_follow_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helicopter_follow_P.PitchTransferFcn_D *
    helicopter_follow_B.PitchCounttorad;

  /* Gain: '<S9>/Gain' */
  helicopter_follow_B.Gain_b = helicopter_follow_P.Gain_Gain_ae * rtb_Backgain;
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' */
    helicopter_follow_B.ElevationCounttorad =
      helicopter_follow_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S6>/Gain' */
    helicopter_follow_B.Gain_e = helicopter_follow_P.Gain_Gain_lv *
      helicopter_follow_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter_follow_B.Sum = helicopter_follow_B.Gain_e +
      helicopter_follow_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helicopter_follow_P.ElevationTransferFcn_C *
    helicopter_follow_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helicopter_follow_P.ElevationTransferFcn_D *
    helicopter_follow_B.ElevationCounttorad;

  /* Gain: '<S7>/Gain' */
  helicopter_follow_B.Gain_dg = helicopter_follow_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  helicopter_follow_B.Gain1[0] = helicopter_follow_P.Gain1_Gain *
    helicopter_follow_B.Sum1;
  helicopter_follow_B.Gain1[1] = helicopter_follow_P.Gain1_Gain *
    helicopter_follow_B.Gain_d;
  helicopter_follow_B.Gain1[2] = helicopter_follow_P.Gain1_Gain *
    helicopter_follow_B.Gain_i;
  helicopter_follow_B.Gain1[3] = helicopter_follow_P.Gain1_Gain *
    helicopter_follow_B.Gain_b;
  helicopter_follow_B.Gain1[4] = helicopter_follow_P.Gain1_Gain *
    helicopter_follow_B.Sum;
  helicopter_follow_B.Gain1[5] = helicopter_follow_P.Gain1_Gain *
    helicopter_follow_B.Gain_dg;
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *)
      helicopter_follow_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helicopter_follow_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_follow_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helicopter_follow_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[100]) {
      currTimeIndex = 99;
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

    helicopter_follow_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Backgain = pDataValues[currTimeIndex];
        } else {
          rtb_Backgain = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Backgain = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 101;
      }
    }
  }

  /* Sum: '<S5>/Sum' incorporates:
   *  Constant: '<S5>/Vd_bias'
   *  Gain: '<S5>/K_pd'
   *  Gain: '<S5>/K_pp'
   *  Sum: '<S5>/Sum2'
   *  Sum: '<S5>/Sum3'
   */
  helicopter_follow_B.Sum_k = ((rtb_Backgain - helicopter_follow_B.Gain1[2]) *
    helicopter_follow_P.K_pp - helicopter_follow_P.K_pd *
    helicopter_follow_B.Gain1[3]) + helicopter_follow_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
  }

  /* Integrator: '<S3>/Integrator'
   *
   * Regarding '<S3>/Integrator':
   *  Limited Integrator
   */
  if (helicopter_follow_X.Integrator_CSTATE >=
      helicopter_follow_P.Integrator_UpperSat ) {
    helicopter_follow_X.Integrator_CSTATE =
      helicopter_follow_P.Integrator_UpperSat;
  } else if (helicopter_follow_X.Integrator_CSTATE <=
             (helicopter_follow_P.Integrator_LowerSat) ) {
    helicopter_follow_X.Integrator_CSTATE =
      (helicopter_follow_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter_follow_X.Integrator_CSTATE;

  /* Sum: '<S3>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   */
  rtb_Derivative = helicopter_follow_P.elevation_ref_Value -
    helicopter_follow_B.Gain1[4];

  /* Sum: '<S3>/Sum2' incorporates:
   *  Constant: '<S3>/Vs_bias'
   *  Gain: '<S3>/K_ed'
   *  Gain: '<S3>/K_ep'
   *  Sum: '<S3>/Sum1'
   */
  helicopter_follow_B.Sum2 = ((helicopter_follow_P.K_ep * rtb_Derivative +
    rtb_Backgain) - helicopter_follow_P.K_ed * helicopter_follow_B.Gain1[5]) +
    helicopter_follow_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helicopter_follow_B.Sum2 - helicopter_follow_B.Sum_k) *
    helicopter_follow_P.Backgain_Gain;

  /* Gain: '<S3>/K_ei' */
  helicopter_follow_B.K_ei = helicopter_follow_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helicopter_follow_DW.TimeStampA >= helicopter_follow_M->Timing.t[0]) &&
      (helicopter_follow_DW.TimeStampB >= helicopter_follow_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Derivative = helicopter_follow_DW.TimeStampA;
    lastU = &helicopter_follow_DW.LastUAtTimeA;
    if (helicopter_follow_DW.TimeStampA < helicopter_follow_DW.TimeStampB) {
      if (helicopter_follow_DW.TimeStampB < helicopter_follow_M->Timing.t[0]) {
        rtb_Derivative = helicopter_follow_DW.TimeStampB;
        lastU = &helicopter_follow_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter_follow_DW.TimeStampA >= helicopter_follow_M->Timing.t[0]) {
        rtb_Derivative = helicopter_follow_DW.TimeStampB;
        lastU = &helicopter_follow_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter_follow_B.PitchCounttorad - *lastU) /
      (helicopter_follow_M->Timing.t[0] - rtb_Derivative);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S10>/Gain' */
  helicopter_follow_B.Gain_l = helicopter_follow_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
  }

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter_follow_P.BackmotorSaturation_UpperSat) {
    helicopter_follow_B.BackmotorSaturation =
      helicopter_follow_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helicopter_follow_P.BackmotorSaturation_LowerSat) {
    helicopter_follow_B.BackmotorSaturation =
      helicopter_follow_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter_follow_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Derivative = (helicopter_follow_B.Sum_k + helicopter_follow_B.Sum2) *
    helicopter_follow_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Derivative > helicopter_follow_P.FrontmotorSaturation_UpperSat) {
    helicopter_follow_B.FrontmotorSaturation =
      helicopter_follow_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Derivative < helicopter_follow_P.FrontmotorSaturation_LowerSat)
  {
    helicopter_follow_B.FrontmotorSaturation =
      helicopter_follow_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter_follow_B.FrontmotorSaturation = rtb_Derivative;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter_follow/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter_follow_DW.HILWriteAnalog_Buffer[0] =
        helicopter_follow_B.FrontmotorSaturation;
      helicopter_follow_DW.HILWriteAnalog_Buffer[1] =
        helicopter_follow_B.BackmotorSaturation;
      result = hil_write_analog(helicopter_follow_DW.HILInitialize_Card,
        helicopter_follow_P.HILWriteAnalog_channels, 2,
        &helicopter_follow_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter_follow_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter_follow_DW.TimeStampA == (rtInf)) {
    helicopter_follow_DW.TimeStampA = helicopter_follow_M->Timing.t[0];
    lastU = &helicopter_follow_DW.LastUAtTimeA;
  } else if (helicopter_follow_DW.TimeStampB == (rtInf)) {
    helicopter_follow_DW.TimeStampB = helicopter_follow_M->Timing.t[0];
    lastU = &helicopter_follow_DW.LastUAtTimeB;
  } else if (helicopter_follow_DW.TimeStampA < helicopter_follow_DW.TimeStampB)
  {
    helicopter_follow_DW.TimeStampA = helicopter_follow_M->Timing.t[0];
    lastU = &helicopter_follow_DW.LastUAtTimeA;
  } else {
    helicopter_follow_DW.TimeStampB = helicopter_follow_M->Timing.t[0];
    lastU = &helicopter_follow_DW.LastUAtTimeB;
  }

  *lastU = helicopter_follow_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter_follow_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter_follow_M->solverInfo);
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
  if (!(++helicopter_follow_M->Timing.clockTick0)) {
    ++helicopter_follow_M->Timing.clockTickH0;
  }

  helicopter_follow_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter_follow_M->solverInfo);

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
    if (!(++helicopter_follow_M->Timing.clockTick1)) {
      ++helicopter_follow_M->Timing.clockTickH1;
    }

    helicopter_follow_M->Timing.t[1] = helicopter_follow_M->Timing.clockTick1 *
      helicopter_follow_M->Timing.stepSize1 +
      helicopter_follow_M->Timing.clockTickH1 *
      helicopter_follow_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter_follow_derivatives(void)
{
  XDot_helicopter_follow_T *_rtXdot;
  _rtXdot = ((XDot_helicopter_follow_T *) helicopter_follow_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_follow_P.TravelTransferFcn_A *
    helicopter_follow_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_follow_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_follow_P.PitchTransferFcn_A *
    helicopter_follow_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_follow_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter_follow_P.ElevationTransferFcn_A *
    helicopter_follow_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter_follow_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S3>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter_follow_X.Integrator_CSTATE <=
            (helicopter_follow_P.Integrator_LowerSat) );
    usat = ( helicopter_follow_X.Integrator_CSTATE >=
            helicopter_follow_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter_follow_B.K_ei > 0)) ||
        (usat && (helicopter_follow_B.K_ei < 0)) ) {
      ((XDot_helicopter_follow_T *) helicopter_follow_M->ModelData.derivs)
        ->Integrator_CSTATE = helicopter_follow_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter_follow_T *) helicopter_follow_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter_follow_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_follow/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter_follow_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (helicopter_follow_DW.HILInitialize_Card,
       "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter_follow_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
      return;
    }

    if ((helicopter_follow_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter_follow_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter_follow_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] =
            helicopter_follow_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter_follow_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] =
            helicopter_follow_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges
        (helicopter_follow_DW.HILInitialize_Card,
         helicopter_follow_P.HILInitialize_analog_input_chan, 8U,
         &helicopter_follow_DW.HILInitialize_AIMinimums[0],
         &helicopter_follow_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_follow_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter_follow_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter_follow_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] =
            helicopter_follow_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter_follow_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] =
            helicopter_follow_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges
        (helicopter_follow_DW.HILInitialize_Card,
         helicopter_follow_P.HILInitialize_analog_output_cha, 8U,
         &helicopter_follow_DW.HILInitialize_AOMinimums[0],
         &helicopter_follow_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_follow_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter_follow_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_follow_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter_follow_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter_follow_DW.HILInitialize_Card,
        helicopter_follow_P.HILInitialize_analog_output_cha, 8U,
        &helicopter_follow_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_follow_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_follow_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter_follow_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter_follow_DW.HILInitialize_Card,
         helicopter_follow_P.HILInitialize_analog_output_cha, 8U,
         &helicopter_follow_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_follow_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter_follow_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter_follow_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter_follow_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter_follow_DW.HILInitialize_Card,
         helicopter_follow_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter_follow_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_follow_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter_follow_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter_follow_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter_follow_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter_follow_DW.HILInitialize_Card,
        helicopter_follow_P.HILInitialize_encoder_channels, 8U,
        &helicopter_follow_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_follow_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter_follow_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_follow_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_follow_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter_follow_DW.HILInitialize_Card,
        helicopter_follow_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter_follow_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter_follow_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues =
          &helicopter_follow_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter_follow_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter_follow_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = helicopter_follow_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter_follow_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = p_HILInitialize_pwm_channels[i1];
            helicopter_follow_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              helicopter_follow_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter_follow_DW.HILInitialize_Card,
          &helicopter_follow_DW.HILInitialize_POSortedChans[0],
          num_duty_cycle_modes,
          &helicopter_follow_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter_follow_DW.HILInitialize_Card,
          &helicopter_follow_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter_follow_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_follow_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            helicopter_follow_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter_follow_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter_follow_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter_follow_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter_follow_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter_follow_DW.HILInitialize_Card,
        helicopter_follow_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *)
        &helicopter_follow_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter_follow_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter_follow_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicopter_follow_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] =
            helicopter_follow_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter_follow_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_follow_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter_follow_DW.HILInitialize_Card,
        helicopter_follow_P.HILInitialize_pwm_channels, 8U,
        &helicopter_follow_DW.HILInitialize_POSortedFreqs[0],
        &helicopter_follow_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_follow_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter_follow_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_follow_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_follow_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter_follow_DW.HILInitialize_Card,
        helicopter_follow_P.HILInitialize_pwm_channels, 8U,
        &helicopter_follow_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_follow_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_follow_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_follow_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter_follow_DW.HILInitialize_Card,
         helicopter_follow_P.HILInitialize_pwm_channels, 8U,
         &helicopter_follow_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter_follow/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (helicopter_follow_DW.HILInitialize_Card,
       helicopter_follow_P.HILReadEncoderTimebase_samples_,
       helicopter_follow_P.HILReadEncoderTimebase_channels, 3,
       &helicopter_follow_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
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
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559753309,
      0.523598775597397, 0.52359877559722157, 0.52359877559699,
      0.523598775596675, 0.52359877559622914, 0.52359877559556334,
      0.5235987755944882, 0.5235987755925201, 0.52359877558795276,
      0.52359877556741929, 0.38856036545868583, 0.1094957602435343,
      -0.11003563544112584, -0.276884522423022, -0.39785602260923381,
      -0.4795618037034845, -0.52358628562450082, -0.52359877490890916,
      -0.523598774954063, -0.52359863675677765, -0.50324239318874853,
      -0.46475678334876425, -0.42054051296355505, -0.37322907148973755,
      -0.32496743350908114, -0.27745317917432505, -0.23198382899978742,
      -0.18950600412492477, -0.15066447991033627, -0.11584962973923739,
      -0.085242140648529391, -0.0588542170351331, -0.036566772640272895,
      -0.018162346496149957, -0.0033536691243817923, 0.0081920446834722258,
      0.016832886062451407, 0.022935851914134192, 0.026863924747256188,
      0.028966268700267188, 0.029571161680182003, 0.028981281303823183,
      0.027470970409997221, 0.02528512286644535, 0.022639348876530574,
      0.019721097611529972, 0.016691430137389291, 0.013687143227252438,
      0.010822940222778845, 0.00819332401767213, 0.0058738464761366512,
      0.0039212915160998525, 0.0023723171285108913, 0.0012401012706759955,
      0.00050879864505998956, 0.00012652119243539912, -2.8712725560718135E-14,
      -1.077019601325731E-14, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter_follow_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_follow_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_follow_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter_follow_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter_follow_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter_follow_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  helicopter_follow_X.Integrator_CSTATE = helicopter_follow_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter_follow_DW.TimeStampA = (rtInf);
  helicopter_follow_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter_follow_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_follow/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter_follow_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter_follow_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter_follow_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter_follow_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_follow_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter_follow_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter_follow_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter_follow_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_follow_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_follow_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter_follow_DW.HILInitialize_Card
                         , helicopter_follow_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter_follow_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter_follow_DW.HILInitialize_AOVoltages[0]
                         , &helicopter_follow_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (helicopter_follow_DW.HILInitialize_Card,
             helicopter_follow_P.HILInitialize_analog_output_cha,
             num_final_analog_outputs,
             &helicopter_follow_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter_follow_DW.HILInitialize_Card,
            helicopter_follow_P.HILInitialize_pwm_channels,
            num_final_pwm_outputs, &helicopter_follow_DW.HILInitialize_POValues
            [0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_follow_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter_follow_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter_follow_DW.HILInitialize_Card);
    hil_close(helicopter_follow_DW.HILInitialize_Card);
    helicopter_follow_DW.HILInitialize_Card = NULL;
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
  helicopter_follow_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter_follow_update();
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
  helicopter_follow_initialize();
}

void MdlTerminate(void)
{
  helicopter_follow_terminate();
}

/* Registration function */
RT_MODEL_helicopter_follow_T *helicopter_follow(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter_follow_P.Integrator_UpperSat = rtInf;
  helicopter_follow_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter_follow_M, 0,
                sizeof(RT_MODEL_helicopter_follow_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter_follow_M->solverInfo,
                          &helicopter_follow_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter_follow_M->solverInfo, &rtmGetTPtr
                (helicopter_follow_M));
    rtsiSetStepSizePtr(&helicopter_follow_M->solverInfo,
                       &helicopter_follow_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter_follow_M->solverInfo,
                 &helicopter_follow_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter_follow_M->solverInfo, (real_T **)
                         &helicopter_follow_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter_follow_M->solverInfo,
      &helicopter_follow_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter_follow_M->solverInfo, (&rtmGetErrorStatus
      (helicopter_follow_M)));
    rtsiSetRTModelPtr(&helicopter_follow_M->solverInfo, helicopter_follow_M);
  }

  rtsiSetSimTimeStep(&helicopter_follow_M->solverInfo, MAJOR_TIME_STEP);
  helicopter_follow_M->ModelData.intgData.f[0] =
    helicopter_follow_M->ModelData.odeF[0];
  helicopter_follow_M->ModelData.contStates = ((real_T *) &helicopter_follow_X);
  rtsiSetSolverData(&helicopter_follow_M->solverInfo, (void *)
                    &helicopter_follow_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter_follow_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter_follow_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter_follow_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter_follow_M->Timing.sampleTimes =
      (&helicopter_follow_M->Timing.sampleTimesArray[0]);
    helicopter_follow_M->Timing.offsetTimes =
      (&helicopter_follow_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter_follow_M->Timing.sampleTimes[0] = (0.0);
    helicopter_follow_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter_follow_M->Timing.offsetTimes[0] = (0.0);
    helicopter_follow_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter_follow_M, &helicopter_follow_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter_follow_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter_follow_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter_follow_M, 25.0);
  helicopter_follow_M->Timing.stepSize0 = 0.002;
  helicopter_follow_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter_follow_M->Sizes.checksums[0] = (3917215120U);
  helicopter_follow_M->Sizes.checksums[1] = (1938816865U);
  helicopter_follow_M->Sizes.checksums[2] = (3435505708U);
  helicopter_follow_M->Sizes.checksums[3] = (1124461689U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter_follow_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter_follow_M->extModeInfo,
      &helicopter_follow_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter_follow_M->extModeInfo,
                        helicopter_follow_M->Sizes.checksums);
    rteiSetTPtr(helicopter_follow_M->extModeInfo, rtmGetTPtr(helicopter_follow_M));
  }

  helicopter_follow_M->solverInfoPtr = (&helicopter_follow_M->solverInfo);
  helicopter_follow_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter_follow_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter_follow_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter_follow_M->ModelData.blockIO = ((void *) &helicopter_follow_B);

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      helicopter_follow_B.Gain1[i] = 0.0;
    }

    helicopter_follow_B.TravelCounttorad = 0.0;
    helicopter_follow_B.Gain = 0.0;
    helicopter_follow_B.Sum1 = 0.0;
    helicopter_follow_B.Gain_d = 0.0;
    helicopter_follow_B.PitchCounttorad = 0.0;
    helicopter_follow_B.Gain_i = 0.0;
    helicopter_follow_B.Gain_b = 0.0;
    helicopter_follow_B.ElevationCounttorad = 0.0;
    helicopter_follow_B.Gain_e = 0.0;
    helicopter_follow_B.Sum = 0.0;
    helicopter_follow_B.Gain_dg = 0.0;
    helicopter_follow_B.Sum_k = 0.0;
    helicopter_follow_B.Sum2 = 0.0;
    helicopter_follow_B.K_ei = 0.0;
    helicopter_follow_B.Gain_l = 0.0;
    helicopter_follow_B.BackmotorSaturation = 0.0;
    helicopter_follow_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter_follow_M->ModelData.defaultParam = ((real_T *)&helicopter_follow_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter_follow_X;
    helicopter_follow_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter_follow_X, 0,
                  sizeof(X_helicopter_follow_T));
  }

  /* states (dwork) */
  helicopter_follow_M->ModelData.dwork = ((void *) &helicopter_follow_DW);
  (void) memset((void *)&helicopter_follow_DW, 0,
                sizeof(DW_helicopter_follow_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_follow_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_follow_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_follow_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_follow_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_follow_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_follow_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_follow_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_follow_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter_follow_DW.TimeStampA = 0.0;
  helicopter_follow_DW.LastUAtTimeA = 0.0;
  helicopter_follow_DW.TimeStampB = 0.0;
  helicopter_follow_DW.LastUAtTimeB = 0.0;
  helicopter_follow_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter_follow_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter_follow_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter_follow_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter_follow_M->Sizes.numY = (0);/* Number of model outputs */
  helicopter_follow_M->Sizes.numU = (0);/* Number of model inputs */
  helicopter_follow_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter_follow_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter_follow_M->Sizes.numBlocks = (56);/* Number of blocks */
  helicopter_follow_M->Sizes.numBlockIO = (18);/* Number of block outputs */
  helicopter_follow_M->Sizes.numBlockPrms = (142);/* Sum of parameter "widths" */
  return helicopter_follow_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
