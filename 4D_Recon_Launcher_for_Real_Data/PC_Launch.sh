
#!/bin/sh -l

export PARAM_INDEX=0
python XT_Main.py --run_setup --run_recon --FBP --ATT --SIM_DATA --num_nodes 1 --PC --Path2Data /blah/ --Path2WhitesDarks /blah/ --Path2Phantom /Users/aditya/Academics/Graduate_Courses/ECE699/Time_Varying_XRay_Tomography/C_Code/Workspace_Argonne/Sim_Datasets/CH_Phantom_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --Path2Mask /Users/aditya/Academics/Graduate_Courses/ECE699/Time_Varying_XRay_Tomography/C_Code/Workspace_Argonne/Sim_Datasets/CH_Mask_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --run_folder ../../Recon_Runs/Att_Recon_Sim_256/ --rot_center 132 --vox_size 0.65 --proj_start 0 --proj_num 1024 --x_width 1024 --recon_x_width 256 --z_start 0 --z_width 4 --recon_z_width 4 --vox_stop_thresh 0.05 --cost_stop_thresh 0.1 --sigma_s 10000 --sigma_t 10000 --K 1 --N_theta 256 --r 1 --multres_xy 3 --multres_z 1 --do_VarEstimate 0 --MaxIter 1000 --msg_string _const_z --min_time_btw_views 0 --rotation_speed 100 --ZingerT 4 --maxHU 40000 --minHU 0           

#python XT_Main.py --run_setup --run_recon --FBP --ATT --SIM_DATA --num_nodes 1 --PC --Path2Data /blah/ --Path2WhitesDarks /blah/ --Path2Phantom /Users/aditya/Academics/Graduate_Courses/ECE699/Time_Varying_XRay_Tomography/C_Code/Workspace_Argonne/Sim_Datasets/CH_Phantom_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --Path2Mask /Users/aditya/Academics/Graduate_Courses/ECE699/Time_Varying_XRay_Tomography/C_Code/Workspace_Argonne/Sim_Datasets/CH_Mask_N_xy_1024_N_z_4_N_p_1024_across_slice_const.bin --run_folder ../../Recon_Runs/Att_Recon_Sim_256/ --rot_center 132 --vox_size 0.65 --proj_start 0 --proj_num 1024 --x_width 1024 --recon_x_width 256 --z_start 0 --z_width 4 --recon_z_width 4 --vox_stop_thresh 0.05 --cost_stop_thresh 0.1 --sigma_s 10000 --sigma_t 10000 --K 16 --N_theta 256 --r 16 --multres_xy 3 --multres_z 1 --do_VarEstimate 0 --MaxIter 1000 --msg_string _const_z --min_time_btw_views 0 --rotation_speed 100 --ZingerT 4 --maxHU 40000 --minHU 0           
