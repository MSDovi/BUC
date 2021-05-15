clear; clc;
tic;
%% Housekeeping
addpath('../aux')
addpath('../aux/svblock')
addpath('../aux/sblock')
addpath('../aux/rsblock')
addpath('../aux/arblock')
addpath('../aux/regblock')
data_raw = readmatrix('../data/GDPC1.csv');
data = 100*log(data_raw(1:270, :));
y = data(:, 2);
T = length(y);

nsim = 10000;
burnin = 1000;
store_theta = zeros(nsim, 5); % [phi, sigc2, sigtau2, tau0]
store_tau = zeros(nsim,T); 
store_h = zeros(nsim,T);
store_g = zeros(nsim,T);

%% Initialise ARBlock
opt_phi.pmean = [1.34 -.7]'; opt_phi.pVar = speye(2);
opt_phi.nu = 3; opt_phi.S = 1*(opt_phi.nu-1);
ar_block = ARBlock(y, zeros(T, 1), zeros(2, 1), opt_phi);

%% Initialise RegBlock
opt_mu.pmean = 2; opt_mu.pVar = 10;
opt_mu.nu = 3; opt_mu.S = 1*(opt_mu.nu - 1);
reg_block = RegBlock(y, ones(T, 1), opt_mu);

%% Initialise State Variables Block
opt_c.ic0 = [0 0]'; opt_c.icV = 10*speye(2);
s_block_c = SBlock(y, zeros(T, 1), opt_phi.pmean, ...
                   reg_block.Var, ar_block.Var, opt_c);
               
opt_tau.ic0 = [y(1)]; opt_tau.icV = 100*speye(1);
rs_block_tau = RSBlock(opt_mu.pmean*ones(T, 1), 1, reg_block.Var, opt_tau);

%% Initialise SVBlocks
opt_sv.nu = 3; opt_sv.S = 0.2; opt_sv.Vs = 0.1;
opt_sv.ic0 = 0; opt_sv.icV = 100;
opt_sv.ncp = true;
sv_block_h = SVBlock(y - opt_mu.pmean, opt_sv);
sv_block_g = SVBlock(y - opt_mu.pmean, opt_sv);

%% Initialise some things defined purely for readability
taum1 = y; % defined only for subsequent readability
Sigma = sparse(1:T,1:T,exp(sv_block_h.svsv));
Omega = sparse(1:T,1:T,exp(sv_block_g.svsv));

%% Gibbs Sampler
for isim = 1:nsim+burnin     
    % Update s_block_c
    s_block_c = s_ex_coeff_var_update(s_block_c, zeros(T, 1), ar_block.coeff, Sigma, Omega);
    s_block_c = s_dep_update(s_block_c, y - reg_block.coeff - taum1);
    s_block_c = s_gibbs_update(s_block_c);
    
    % Update rs_block_tau
    rs_block_tau = rs_ex_coeff_var_update(rs_block_tau, reg_block.coeff*ones(T, 1), 1, Sigma);
    rs_block_tau = rs_s_update(rs_block_tau, y - s_block_c.s); % Force change tau because it is a residual
    rs_block_tau = rs_gibbs_update(rs_block_tau);
    % unpack variables for easier reading
    tau = y - s_block_c.s;
    taum1 = [rs_block_tau.ic; tau(1:end-1)];
    
    % Update reg_bock
    reg_block = reg_dep_indie_update(reg_block, tau - taum1, ones(T, 1));
    reg_block = reg_var_update(reg_block, Sigma);
    reg_block = reg_gibbs_update(reg_block);
    
    % Update ar_block
    ar_block = ar_dep_ic_update(ar_block, s_block_c.s, s_block_c.ic);
    ar_block = ar_ex_var_update(ar_block, zeros(T, 1), Omega);
    ar_block = ar_gibbs_update(ar_block);

    % Update sv_block_h and sv_block_g
    sv_block_h = sv_dep_update(sv_block_h, rs_block_tau.res);
    sv_block_h = sv_gibbs_update(sv_block_h);    
    sv_block_g = sv_dep_update(sv_block_g, s_block_c.res);
    sv_block_g = sv_gibbs_update(sv_block_g);
    % unpack variables for easier reading
    Sigma = sparse(1:T,1:T,exp(sv_block_h.svsv));
    Omega = sparse(1:T,1:T,exp(sv_block_g.svsv));
        
    if (mod(isim, 1000) ==0)
        disp([num2str(isim) ' loops... ']);
    end
    
    if isim > burnin
        i = isim-burnin;
        store_tau(i,:) = tau';
        store_theta(i,:) = [ar_block.coeff' ar_block.Var(1,1) reg_block.Var(1,1) rs_block_tau.ic'];
        store_h(i,:) = sv_block_h.svsv'; 
        store_g(i,:) = sv_block_g.svsv'; 
    end    
end
       
tau_hat = mean(store_tau)';
theta_hat = mean(store_theta)'
theta_CI = quantile(store_theta,[.025 .975])
h_hat = mean(exp(store_h/2))'; 
g_hat = mean(exp(store_g/2))'; 

%% plot of graphs
%tt = (1980.25:2019.5)';
tt = (1:158)';
figure;
hold on
    plot(y-tau_hat,'linewidth',1);
    plot(zeros(T,1),'--k','linewidth',1);
hold off
title('Output gap estimates');
box off;
set(gcf,'Position',[100 100 800 300]);

figure;
plot([tau_hat, y]); box off;
set(gcf,'Position',[100 100 800 300]);
title('Trend output growth estimates');

figure;
plot([h_hat g_hat]); box off;
set(gcf,'Position',[100 100 800 300]);