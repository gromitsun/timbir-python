#!/bin/bash

START_TIME=$(date +%s)

./XT_Engine --p 1.2 --sigma_s 500000 --sigma_t 50000 --c_s 0.0000001 --c_t 0.0000001 --delta_xy 4 --delta_z 4 --length_r 1000 --length_t 62.5 --voxel_thresh 5 --cost_thresh 0.1 --iter 1000 --rotation_center 128 --alpha 1.5 --sinobin 2 --initICD 0 --Rtime0 0 --Rtime_delta 64 --Rtime_num 16 --num_projections 1024 --N_r 256 --N_t 16 --detector_slice_begin 2 --detector_slice_end 13 --phantom_N_xy 512 --phantom_N_z 16 --num_threads 8 --time_reg --writeTiff
#./XT_Engine --p 1.2 --sigma_s 500000 --sigma_t 50000 --c_s 0.0000001 --c_t 0.0000001 --delta_xy 2 --delta_z 2 --length_r 1000 --length_t 62.5 --voxel_thresh 5 --cost_thresh 0.1 --iter 1000 --rotation_center 128 --alpha 1.5 --sinobin 2 --initICD 0 --Rtime0 0 --Rtime_delta 64 --Rtime_num 16 --num_projections 1024 --N_r 256 --N_t 16 --detector_slice_begin 2 --detector_slice_end 13 --phantom_N_xy 512 --phantom_N_z 16 --num_threads 8 --time_reg --writeTiff
#./XT_Engine --p 1.2 --sigma_s 500000 --sigma_t 50000 --c_s 0.0000001 --c_t 0.0000001 --delta_xy 1 --delta_z 1 --length_r 1000 --length_t 62.5 --voxel_thresh 5 --cost_thresh 0.1 --iter 1000 --rotation_center 128 --alpha 1.5 --sinobin 2 --initICD 0 --Rtime0 0 --Rtime_delta 64 --Rtime_num 16 --num_projections 1024 --N_r 256 --N_t 16 --detector_slice_begin 2 --detector_slice_end 13 --phantom_N_xy 512 --phantom_N_z 16 --num_threads 8 --time_reg --writeTiff

END_TIME=$(date +%s)
DIFF_TIME=$(( $((END_TIME-START_TIME))/60 ))
echo "run: Total execution time of script is $DIFF_TIME minutes"



