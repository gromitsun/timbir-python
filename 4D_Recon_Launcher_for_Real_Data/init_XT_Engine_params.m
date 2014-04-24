
function [param] = init_XT_Engine_params ()

param.views.delta_tau = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.num_threads = 8;
param.reconstruct = 1;

param.voxel_thresh = [5 5 5 5];
param.cost_thresh = [0.1 0.1 0.1 0.1];
param.delta_x = [8 4 2 1];
param.initICD = [0 2 2 2];
param.writeTiff = [0 1 1 1];

param.c_s = 10^-6;
param.c_t = 10^-6;

param.N_r = 512;
param.N_p = 1280;

% param.r = [1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
% param.views.K = [256 256 256 256 1 1 1 1 64 64 64 64 16 16 16 16 4 4 4 4];
% param.views.L = [1 1 1 1 256 256 256 256 4 4 4 4 16 16 16 16 64 64 64 64];

param.r = [16 1];
param.views.K = [16 1];
param.views.L = 2048./param.views.K;

% param.sigma_s = [5*10^5];
% param.sigma_t = [5*10^3];

sigma_s = [10^5];
sigma_t = [10^4];
% [sig_s,sig_t] = meshgrid(sigma_s, sigma_t);
for i=1:numel(param.r)    
    param.regp(i).sigma_s = sigma_s;
    param.regp(i).sigma_t = sigma_t;
end

%param.comparison_times = (128:128*(16-1)-1)*param.views.delta_tau;

param.comparison_times = (128:128*(10-1)-1)*param.views.delta_tau;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%param.gen_projections_method = 'real_data';
param.gen_projections_method = 'phantom';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.length_r = 1500;
param.views.time0 = 0;

param.projection.radius = param.length_r/4;
param.projection.time_start = param.views.time0;
param.projection.time_stop = param.comparison_times(numel(param.comparison_times)/2);

param.add_noise_to_projection = 1;
param.views.N_theta = (param.views.K).*(param.views.L);
param.copy_source = 0;
param.copy_executables = 0;
param.use_tifflibrary = 0;

param.p = 1.2;
param.iter = 1000;
param.rotation_center = param.N_r/2;
param.alpha = 1.5;
param.time_reg = 1;
param.Rtime0 = 0;
param.Rtime_num = floor(param.r.*(param.N_p./param.views.N_theta));
param.Rtime_delta = (param.N_p*param.views.delta_tau)./param.Rtime_num;

param.multi_res_stages = numel(param.delta_x);

param.N_x = param.delta_x(end)*param.N_r;
param.N_y = param.N_x;

if(numel(param.r) ~= numel(param.views.K) || numel(param.r) ~= numel(param.views.L))
    display('ERROR: array length dont match');
    return;
end

param.maxHU = 43000;
param.minHU = 7000;
