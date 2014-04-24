function [views, buf1, buf2] = gen_interlaced_views_0_to_inf (L, K, N_p)

k = log2(K);
N_theta = L*K;
delta_theta = 180/N_theta;

views = zeros(N_p, 1);
buf1 = zeros(N_p, 1);
buf2 = zeros(N_p, 1);
for i=0:N_p-1
   buf1(i+1) = i*K;
   buf2(i+1) = bitreversed_decimal(mod(floor(i/L), K), k);
   views(i+1) = buf1(i+1) + buf2(i+1); 
   views(i+1) = views(i+1)*delta_theta;
end


