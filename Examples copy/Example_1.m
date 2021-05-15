clear; clc;
tic;
%% Housekeeping
addpath('../aux');
addpath('../aux/svblock')
addpath('../aux/regblock')
y = load('../data/AUDUSD.csv');
T = length(y);

nsim = 20000;
burnin = 5000;
store_theta = zeros(nsim,3); 
store_h = zeros(nsim,T);

%% Initalise the Regblock for mu
opt_mu.pmean = 0; opt_mu.pVar = 100;
opt_mu.nu = 3; opt_mu.S = 0.2;

obs_reg_block = RegBlock(y, ones(T, 1), opt_mu);

%% Initialise the SVBlock for h
opt_sv.nu = 3; opt_sv.S = 0.2; opt_sv.Vs = 0.1;
opt_sv.ic0 = 0; opt_sv.icV = 100;
opt_sv.ncp = false;

sv_block = SVBlock(y, opt_sv);

%% Gibbs Sampler
for isim = 1:nsim + burnin
    
    % Update reg_block
    obs_reg_block = reg_var_update(obs_reg_block, ...
                                   sparse(1:T,1:T,exp(sv_block.svsv)));
    obs_reg_block = reg_gibbs_update(obs_reg_block);    
    
    % Update sv_block}
    sv_block = sv_dep_update(sv_block, obs_reg_block.res);
    sv_block = sv_gibbs_update(sv_block);

    if (mod(isim, 5000) == 0)
        disp([num2str(isim) ' loops... ']);
    end    
    
    if isim > burnin
        isave = isim - burnin;
        store_h(isave, :) = sv_block.svsv';
        store_theta(isave,:) = [obs_reg_block.coeff sv_block.ic ... 
                                                          sv_block.var];
    end    
end

%% Post-Simulation Manipulation
theta_hat = mean(store_theta);
theta_CI = quantile(store_theta,[.025 .975]);
h_hat = mean(exp(store_h/2))'; 

%% Graphing
tt = linspace(2005,2013,T)';
figure;
plot(tt,h_hat); box off; xlim([2005 2013]);
set(gcf,'Position',[100 100 800 300]);
toc;