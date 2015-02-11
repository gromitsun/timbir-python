#!/bin/sh -l

#python XT_Extract_Slices_4m_Dataset.py --Path2Source $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/APS14_60_AlSi.hdf --Path2Dest $RCAC_SCRATCH/Argonne_Datasets/APS_Beamtime_042014/AlSi_data_2.hdf --z_start_extract 300 --z_num_extract 4 --contain_darks

python XT_Extract_Slices_4m_Dataset.py --Path2Source /projects/b1003/Ashwin/APS_Aug2014/Datasets/SlowCool_FT_606to591_2C_CR_2ms_ET.hdf --Path2Dest /projects/b1003/Ashwin/APS_Aug2014/Datasets/SlowCool_FT_606to591_2C_CR_2ms_ET_extract.hdf --z_start_extract 300 --z_num_extract 50 --contain_whites --contain_darks

#python XT_SubViews_ExtractSlices_4m_Dataset.py --Path2Source /projects/p20636/SlowCool_FT_612to602.hdf --Path2Dest /projects/p20636/SlowCool_FT_612to602_sub4.hdf --z_start_extract 500 --z_num_extract 8 --view_subsmpl 4 --contain_whites --contain_darks

