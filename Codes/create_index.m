clc;
clear all;

mu_v = [0.01:0.005:0.05];
alpha_v = [1.0:0.05:1.8];
%beta_v = [0.03:0.01:0.09 0.1:0.05:0.25 0.3:0.1:0.5];

n = 0;
for i = 1:numel(mu_v)
    for j = 1:numel(alpha_v)
       % for k = 1:numel(beta_v)
        mu = mu_v(i);
        a = alpha_v(j);
        %beta = beta_v(k);
        n = n+1;
        parameter(n,:) = [n,mu,a]; 
        %end        
    end
end
outdir = pwd;
outpath = fullfile(outdir,'parameter.mat');
save(outpath,'parameter')
