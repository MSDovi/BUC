clear; clc;
tic;
%% Housekeeping
addpath('../aux');
addpath('../aux/svblock')
addpath('../aux/regblock')
data_raw = load('../data/USCPI_2015Q4.csv');
data = 400*log(data_raw(2:end)./data_raw(1:end-1));
y = data;
T = length(y);

nsim = 20000;
burnin = 1000;
valh = 0; valg = 0;

store_theta = zeros(nsim,5); % [omegah omegag h0 g0 tau0]
store_tau = zeros(nsim,T); 
store_h = zeros(nsim,T);
store_g = zeros(nsim,T);
store_lpden = zeros(nsim,3); % log posterior densities of 
                             % omegah = valh; omegag = valg; 
                             % omegah = omegag

%% Initialise the SVBlock for h, g                            
opt_sv.nu = 3; opt_sv.S = 0.2; opt_sv.Vs = 0.1;
opt_sv.ic0 = 0; opt_sv.icV = 10;
opt_sv.ncp = true;

sv_block_h = SVBlock(y, opt_sv);
sv_block_g = SVBlock(y, opt_sv);

%% Initialise the SBlock for tau
opt_tau.ic0 = [0]'; opt_tau.icV = [10];

s_block_tau = SBlock(y, zeros(T, 1), 1, ... 
                              sparse(1:T,1:T,exp(sv_block_h.svsv)), ... 
                              sparse(1:T,1:T,exp(sv_block_g.svsv)), ...
                              opt_sv);

%% Gibbs Sampler
for isim = 1:nsim+burnin
    
    % Update s_block_tau
    s_block_tau = s_ex_coeff_var_update(s_block_tau, zeros(T, 1), 1, ...
                                sparse(1:T,1:T,exp(sv_block_h.svsv)), ...
                                sparse(1:T,1:T,exp(sv_block_g.svsv)));
    s_block_tau = s_gibbs_update(s_block_tau);
    
    % Update sv_block_h and sv_block_g
    sv_block_h = sv_dep_update(sv_block_h, y - s_block_tau.s);
    sv_block_h = sv_gibbs_update(sv_block_h);    
    sv_block_g = sv_dep_update(sv_block_g, ... 
                    s_block_tau.s-[s_block_tau.ic;s_block_tau.s(1:end-1)]);
    sv_block_g = sv_gibbs_update(sv_block_g);

    if (mod(isim, 5000) == 0)
        disp([num2str(isim) ' loops... ']);
    end     
    
    if isim > burnin
        isave = isim - burnin;
        store_tau(isave,:) = s_block_tau.s';
        store_h(isave,:) = sv_block_h.svsv'; 
        store_g(isave,:) = sv_block_g.svsv'; 
        store_theta(isave,:) = [sv_block_h.sd sv_block_g.sd ... 
                               sv_block_h.ic sv_block_g.ic s_block_tau.ic]; 

        lh0 = -.5*log(2*pi*sv_block_h.Dsd) - ... 
                              .5*(sv_block_h.sd_hat-valh)^2/sv_block_h.Dsd;
        lg0 = -.5*log(2*pi*sv_block_g.Dsd) - ... 
                              .5*(sv_block_g.sd_hat-valg)^2/sv_block_g.Dsd;
        lhg0 = lh0 + lg0;
        store_lpden(isave,:) = [lh0 lg0 lhg0];
       
    end    
end

%% Post-Simulation Manipulation

theta_hat = mean(store_theta)';
tau_hat = mean(store_tau)';
h_hat = mean(exp(store_h/2))'; 
g_hat = mean(exp(store_g/2))'; 

maxlpden = max(store_lpden);
lpostden = log(mean(exp(store_lpden-repmat(maxlpden,nsim,1)))) + maxlpden;
lpriden = [log(normpdf(valh,0,sqrt(sv_block_h.opt.Vs))) ... 
           log(normpdf(valg,0,sqrt(sv_block_g.opt.Vs)))];
lpriden(3) = sum(lpriden(1:2));
lBF = lpriden-lpostden;

%% Graphing
tt = (1947.25:.25:2015.75)';
figure;
plot(tt,[tau_hat y]); xlim([1947 2016]); box off;
set(gcf,'Position',[100 100 800 300]);

figure;
plot(tt,[h_hat g_hat]); xlim([1947 2016]); box off;
set(gcf,'Position',[100 100 800 300]);

figure;
subplot(1,2,1); hist(store_theta(:,1),50); box off;
subplot(1,2,2); hist(store_theta(:,2),50); box off;
set(gcf,'Position',[100 100 800 300]);
toc;