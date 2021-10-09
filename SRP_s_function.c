#define S_FUNCTION_NAME  SRP_s_function
#define S_FUNCTION_LEVEL 2
#include "simstruc.h"

#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "SRP_acc.h"

static void mdlInitializeSizes(SimStruct *S)
{
    // S-function parameters
    // none
   //  ssSetNumSFcnParams(S, 2);
   //  if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
   //     return; /* Parameter mismatch reported by the Simulink engine*/
   // }
    // ******************************************************************
    // * Input ports:                                                   *
    // ******************************************************************
    // * I. Satellite position:                                            *
    // *     1-3 Satellite position                                     *
    // ******************************************************************
    // Set number of input ports to 1
    if (!ssSetNumInputPorts(S, 1)) {
        return;
    }
    // Set up input Ports
    // PORT 1: SATELLITE POSITION
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetInputPortDirectFeedThrough(S, 0, 1); // Port will be used in method mldOutputs.(Direct feedthrough)
    ssSetInputPortRequiredContiguous(S, 0, 1); /* Returns true if the signal elements entering the specified port must occupy
                                                   contiguous areas of memory. If the elements are contiguous, a method can
                                                    access the elements of the signal simply by incrementing the signal pointer
                                                    returned by ssGetInputPortSignal. */
   
   
    // Set up output Ports  
    if (!ssSetNumOutputPorts(S,1)) {
        return;
    }
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);

    // Two outputs with tone entry each: C_D and C_L
    ssSetOutputPortWidth(S, 0, DYNAMICALLY_SIZED);


    
    ssSetNumSampleTimes(S, 1);
    
    /* Initialize DWork vectors to save or hold parameters and variables form previous steps */
   ssSetNumDWork(S, 1);
   ssSetDWorkWidth(S, 0, -1);
   ssSetDWorkDataType(S, 0, SS_DOUBLE); 
   
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, -1);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
}

// Provide input and output port dimension info (otherwise dynmaic arrays do not work)
#define MDL_SET_INPUT_PORT_DIMENSION_INFO
#if defined(MDL_SET_INPUT_PORT_DIMENSION_INFO) && defined(MATLAB_MEX_FILE)

/* Function: mdlSetInputPortDimensionInfo =================================
 * Abstract:
 *    This method is called with the candidate dimensions for an input port
 *    with unknown dimensions. If the proposed dimensions are acceptable, the
 *    method should go ahead and set the actual port dimensions.
 *    If they are unacceptable an error should be generated via
 *    ssSetErrorStatus.
 *    Note that any other input or output ports whose dimensions are
 *    implicitly defined by virtue of knowing the dimensions of the given
 *    port can also have their dimensions set.
 *
 */
static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T portIndex, const
  DimsInfo_T *dimsInfo)
{
/* /* Set input port dimension */
if(!ssSetInputPortDimensionInfo(S, portIndex, dimsInfo)) return;


}
#endif


#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO
#if defined(MDL_SET_OUTPUT_PORT_DIMENSION_INFO) && defined(MATLAB_MEX_FILE)
/* Function: mdlSetOutputPortDimensionInfo ================================
 * Abstract:
 *    This method is called with the candidate dimensions for an output port
 *    with unknown dimensions. If the proposed dimensions are acceptable, the
 *    method should go ahead and set the actual port dimensions.
 *    If they are unacceptable an error should be generated via
 *    ssSetErrorStatus.
 *    Note that any other input or output ports whose dimensions are
 *    implicitly defined by virtue of knowing the dimensions of the given
 *    port can also have their dimensions set.
 *
 */
static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T portIndex, const DimsInfo_T *dimsInfo)
{
/* Set output port dimension */
if(!ssSetOutputPortDimensionInfo(S, portIndex, dimsInfo)) return;


}
#endif

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    Initialize both continuous states to zero
 */
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *x   = (real_T*) ssGetDWork(S,0);


  
}
#define MDL_OUTPUT
static void mdlOutputs(SimStruct *S, int_T tid)
{
    
    int vec_size, n_sat,i;
    double param_arr_in[3], param_arr_out[3];
  
    // Get parameters from inputs
  real_T *sat_position = (real_T *) ssGetInputPortRealSignal(S, 0);
  
  // Get the size of input vector
  
  int *vec_dim  = ssGetInputPortDimensions(S,0);
  // printf("vec: %d \n",vec_dim[0]);
  
  // number of satellites  
  n_sat = vec_dim[0]/3;  

  // Define output variable  
    real_T *new_acc_SRP = (real_T *) ssGetOutputPortRealSignal(S, 0);


    //double output[3];
    for (i=1; i<=n_sat; i++)
    {
        param_arr_in[0] = sat_position[3*i-3];
        param_arr_in[1] = sat_position[3*i-2];
        param_arr_in[2] = sat_position[3*i-1];
        
        // Call c-fucntion
        
        SRP_acc(param_arr_in, param_arr_out);
        
        // Return calculated acceleration vector
        
        new_acc_SRP[3*i-3] = param_arr_out[0];
        new_acc_SRP[3*i-2] = param_arr_out[1];
        new_acc_SRP[3*i-1] = param_arr_out[2];
       
    }
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{

}

static void mdlTerminate(SimStruct *S)
{
 
}
#ifdef MATLAB_MEX_FILE    /* Is this file being compiled as a 
                             MEX-file? */
#include "simulink.c"     /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"      /* Code generation registration 
                             function */
#endif