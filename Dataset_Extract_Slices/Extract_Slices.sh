#!/bin/sh -l

python XT_Extract_Slices_4m_Dataset.py --Path2Source $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_60_AlSi.hdf --Path2Dest $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/AlSi_data_2.hdf --z_start_extract 300 --z_num_extract 4 --contain_darks

