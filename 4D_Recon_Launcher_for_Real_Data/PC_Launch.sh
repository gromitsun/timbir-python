
#!/bin/sh -l

export PARAM_INDEX=1
#python XT_Main.py --setup_launch_folder --run_reconstruction --MBIR_ATT_SIM --num_MPI_process 1 --PC
python XT_Main.py --setup_launch_folder --run_reconstruction --MBIR_PHCON_REAL --num_MPI_process 1 --PC
