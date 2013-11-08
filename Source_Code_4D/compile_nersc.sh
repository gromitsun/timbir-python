cc -Wall -ansi -fopenmp -DEXPECTED_COUNTS_FOR_PHANTOM_DATA="29473" -DPHANTOM_FILENAME="\"phantom_shrinking_sphere\"" -DDATA_HDF_FILENAME="\"data.hdf5\"" -DWHITEDARK_HDF_FILENAME="\"data.hdf5\"" -DPROJECTION_HDF_START="0" -DBH_QUAD_COEF="0.5" -DPOSITIVITY_CONSTRAINT -DNO_COST_CALCULATE -DEXPECTED_COUNTS_FOR_PHANTOM_DATA="29473" -DPHANTOM_FILENAME="\"phantom_shrinking_sphere\"" -o XT_Engine XT_Engine.c XT_ICD_update.c XT_Init.c XT_genSinogram.c XT_AMatrix.c XT_Profile.c XT_NHICD.c allocate.c randlib.c tiff.c XT_IOMisc.c XT_ImageProc.c XT_MPI.c XT_HDFIO.c -lm 


