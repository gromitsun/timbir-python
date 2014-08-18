clear; close all;

K = 8;
r = 8;
N_theta = 256;
N_p = 256*4;
N_r = 256;
N_t = 4;
rotation_center = 131;
node_num = 1;
length_r = 0.65*1024;
path2result = '../../Recon_Runs/Att_Recon_Sim_256/XT_Result_Repository/FBP_r_8_K_8_N_theta_256_N_p_1024/';
path2launch = '../../Recon_Runs/Att_Recon_Sim_256/XT_run/run_sigs_300000.0_sigt_3000.0_r_8_K_8_N_theta_256_N_p_1024_sig2_10_th0.1/';
L = N_theta/K;
[views, buf1, buf2] = gen_interlaced_views_0_to_inf (L, K, N_p);

Rtime_num = r*N_p/N_theta;
Rtime_delta = N_p/Rtime_num;

ratio = N_t/node_num;
centered_proj_sz = 2*rotation_center + 1;
x_idx = 1:min(N_r,centered_proj_sz);
y_idx = (1:N_r) + floor((centered_proj_sz-N_r)/2);
for rank = 1:node_num
    fid = fopen([path2launch,'projection_n', num2str(rank-1), '.bin'], 'r');
    proj_temp = fread(fid, N_p*N_r*ratio,'double');
    fclose(fid);
    projections = permute(reshape(proj_temp,[ratio,N_r,N_p]),[2,3,1]); 
    centered_projections = zeros(centered_proj_sz, N_p, ratio); 
    centered_projections(1:min(centered_proj_sz,N_r), :, :) = projections(1:min(centered_proj_sz,N_r), :, :);

    recon = zeros(N_r, N_r, ratio); 
    for i = 1:Rtime_num
        for j = 1:ratio
            temp = iradon(centered_projections(:,(i-1)*N_theta/r+1:i*N_theta/r,j), views((i-1)*N_theta/r+1:i*N_theta/r), 'cubic', 'Hamming', 1, centered_proj_sz);
            recon(:,1:min(N_r,centered_proj_sz),j) = N_r/length_r*temp(y_idx,1:min(N_r,centered_proj_sz));
        end
        recon(recon < 0) = 0;
        fid = fopen([path2launch, 'object_n', num2str(rank-1), '_time_', num2str(i-1), '.bin'], 'wb');
        fwrite(fid, permute(recon,[2,1,3]), 'double');
        fclose(fid);    
    end
end
