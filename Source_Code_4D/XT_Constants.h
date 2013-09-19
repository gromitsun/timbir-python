
#ifndef XT_CONSTANTS_H
#define XT_CONSTANTS_H

/*#define NO_COST_CALCULATE*/
/*#define DEBUG_EN*/
#define ZERO_SKIPPING
#define INIT_SPARSE_ANGLES_FROM_FILE

/*#define FOLDER_PATH "/home/aditya/Academics/ECE699/"*/
#define FOLDER_PATH "/Users/aditya/Academics/Graduate_Courses/ECE699/Time_Varying_XRay_Tomography/Data/"
#define FBP_INNER_FOLDER "/rec_8bit/"
#define PROJECTION_INNER_FOLDER "/tif/"
#define FILE_FIRSTNAME "26B_"
#define AVG_UPDATE_FILENAME "avg_update"

#define TIME_BEGIN 1
#define TIME_END 104
#define ANGLE_BEGIN 5
#define ANGLE_END 724
#define ROW_BEGIN 512
#define ROW_END 512
#define COL_BEGIN 0
#define COL_END 1023

#define TIME_DIGITS 4
#define ANGLE_DIGITS 4
#define REPEAT_DIGITS "0000"

/* i.e., the number of digits in projection angle part of filename and number of 0 in REPEAT_DIGITS should be greater than max(TIME_DIGITS,ANGLE_DIGITS) */

#define PROJECTION_FILENAME "projection"
#define WEIGHT_MATRIX_FILENAME "weight"
#define BRIGHT_FIELD_FILENAME "bright"
#define OBJECT_FILENAME "object"
#define PROJ_OFFSET_FILENAME "proj_offset"

#ifndef PHANTOM_FILENAME
  #define PHANTOM_FILENAME "phantom_shrinking_sphere"
#endif

#define PHANTOM_FOLDER "/Users/aditya/Academics/Graduate_Courses/ECE699/Time_Varying_XRay_Tomography/C_Code/Matlab_Launcher_Continuous_Views/Phantoms/"
#ifndef SPARSE_ANGLES_LIST_FILE
  #define SPARSE_ANGLES_LIST_FILE "view_info.txt"
#endif

#define OBJECT_INIT_VAL 0

/*#define AIR_MASS_ATT_COEFF 0.5658*/ /*in cm^2/g*/
/*#define WATER_MASS_ATT_COEFF 0.5926*/
#define HFIELD_UNIT_CONV_CONST 0.0001
#define AIR_MASS_ATT_COEFF 0.496372335005353 /*in cm^2/g. computed using cubic interpolation*/
#define WATER_MASS_ATT_COEFF 0.521225397034623

#define WATER_DENSITY 1.0 /*in g/cm^3*/
#define AIR_DENSITY 0.001205
#define HOUNSFIELD_WATER_MAP 1000
#define HOUNSFIELD_AIR_MAP 0
#ifndef HOUNSFIELD_MAX
  #define HOUNSFIELD_MAX 60000
#endif
#ifndef HOUNSFIELD_MIN
  #define HOUNSFIELD_MIN 10000
#endif
#define MRF_Q 2.0

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2         1.57079632679489661923132169163975144   /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4         0.785398163397448309615660845819875721  /* pi/4 */
#endif

#define PROFILE_RESOLUTION 1536
#define BEAM_RESOLUTION 512
#define DETECTOR_RESPONSE_BINS 64
#define QGGMRF_ITER 1
#define TRUE_ANGULAR_STEP_SIZE 0.25

#define NHOOD_Z_MAXDIM 3
#define NHOOD_Y_MAXDIM 3
#define NHOOD_X_MAXDIM 3
#define NHOOD_TIME_MAXDIM 3

#define BRIGHT_FIELD_INDEX_1 3
#define BRIGHT_FIELD_INDEX_2 4

#define PRE_FLAT_DARK_AQ_TIME 1.0
#define MV_SAMPLE_IN_BEAM_TIME 20.0
#define MV_SAMPLE_OUT_BEAM_TIME 20.0
#define POST_FLAT_DARK_AQ_TIME 0.5
#define WAIT_TIME_BTW_CYCLES 11.5
#define PROJECTION_AQ_TIME 180.0
#define BH_QUAD_COEF 0.5

#ifndef EXPECTED_COUNTS_FOR_PHANTOM_DATA
  #define EXPECTED_COUNTS_FOR_PHANTOM_DATA 29473
#endif


#endif /*#ifndef XT_CONSTANTS_H*/
