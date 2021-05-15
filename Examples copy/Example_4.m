%% Housekeeping
clear; clc;
tic;
addpath('../aux')
addpath('../aux/svblock')
addpath('../aux/sblock')
addpath('../aux/arblock')
addpath('../aux/dfblock')
data = load('../data/freddata_Q.csv');
data_mean = mean(data)';
data_std = std(data)';
idx = ~isnan(data_mean); % extract variables that don't have missing data
n = sum(idx);
T = size(data,1);
Y = (data(:,idx) - repmat(data_mean(idx)',T,1))./repmat(data_std(idx)',T,1); % standardised the data
y = reshape(Y',T*n,1);
 
nsim = 1000;
burnin = 1000;
k = 1; % # of factors
 
store_F = zeros(nsim,T,k);
store_Lambda = zeros(nsim,n,k);
store_g = zeros(nsim,T);
store_h = zeros(nsim, T, n);
 
%% Initialise DFBlock
nLambda = (2*n-k-1)*k/2;
opt.ic0 = zeros(k, 1); opt.icV = 0.000000001*eye(k);
opt.ml = zeros(nLambda, 1); opt.Vl = 100 * speye(nLambda);
opt.nu = 5; opt.S = 1*(opt.nu-1);
opt.T = T;
Phi = 0.6;
Lambda_init = zeros(n,k);
Lambda_init(1:n+1:end) = 1; % set diagonal elements to 1

df_block = DFBlock(Y, zeros(T*k, 1), Lambda_init, Phi, sparse(diag(1)), opt);
 
%% Initialise ARBlock
opt_phi.pmean = [0.5]'; opt_phi.pVar = speye(1);
opt_phi.nu = 5; opt_phi.S = 1*(opt_phi.nu-1);
 
ar_block_f = ARBlock(df_block.f, zeros(T, 1), df_block.ic, opt_phi);
 
%% Initialise SVBlocks
opt_sv.nu = 3; opt_sv.S = 0.2; opt_sv.Vs = 0.1;
opt_sv.ic0 = 0; opt_sv.icV = 100;
opt_sv.ncp = true;
 
sv_block_g = SVBlock(df_block.f - df_block.exvar, opt_sv);
 
for i = 1:n
    obs_sv_blocks(i) = SVBlock(df_block.Res_obs(:, i), opt_sv);
end

%% Initialise some things defined purely for readability
flat_Sigma = zeros(T*n, 1);

%% Gibbs Sampler
disp('Starting MCMC for the dynamic factor model.... ');
for isim = 1:nsim + burnin
    isim
    %Unpack some things for easier reading
    %Omega = ar_block_f.Var; use this if you want the homo
    %Sigma = df_block.obsVar;
    Omega = sparse(1:T,1:T,exp(sv_block_g.svsv)); % T x T
    for i = 1:n
       slct = i:n:((T-1)*n+i);
       flat_Sigma(slct) = exp(obs_sv_blocks(i).svsv);
    end
    Sigma = sparse(1:(T*n), 1:(T*n), flat_Sigma); % Tn x Tn

    Phi = ar_block_f.coeff;
        
    % sample f and factor loadings
    df_block = df_ex_coeff_var_update(df_block, zeros(T*k, 1), Phi, Sigma, Omega);
    df_block = df_gibbs_update(df_block);
    % Unpack some things for easier reading
    F = df_block.F; F_0 = df_block.F0;
    Lambda = df_block.loadings_s;
 
    % sample Phi and Omega from ar_block
    ar_block_f = ar_dep_ic_update(ar_block_f, df_block.f, df_block.ic);
    ar_block_f = ar_ex_var_update(ar_block_f, zeros(T, 1), Omega);
    ar_block_f = ar_gibbs_update(ar_block_f);
 
    % SV Block
    for i = 1:n
        obs_sv_blocks(i) = sv_dep_update(obs_sv_blocks(i), df_block.Res_obs(:, i));
        obs_sv_blocks(i) = sv_gibbs_update(obs_sv_blocks(i));
    end
    sv_block_g = sv_dep_update(sv_block_g, df_block.Res);
    sv_block_g = sv_gibbs_update(sv_block_g);
    
    if isim > burnin
        isave = isim - burnin;
        store_F(isave,:,:) = F;
        store_Lambda(isave,:,:) = Lambda;        
        store_g(isave,:) = sv_block_g.svsv'; 
        for i = 1:n
            store_h(isave, :, i) = obs_sv_blocks(i).svsv';
        end
    end
    
     if (mod(isim,2000) == 0)
        disp([num2str(isim) ' loops... '])
    end
    
end
 
F_hat = squeeze(mean(store_F));
Lambda_hat = squeeze(mean(store_Lambda));
g_hat = mean(exp(store_g/2))'; 
h_hat = squeeze(mean(exp(store_h)));

%% plot of graphs
dates = linspace(1959.75,2015.75,T)';
figure;
hold on
    plot(dates,F_hat); 
    hline = refline(0,0);
    hline.Color = 'k';    
hold off
box off; 
xlim([dates(1) dates(end)]);  
set(gcf,'Position',[100 100 800 300]);
 
figure;
plot(dates,g_hat); xlim([dates(1) dates(end)]);  box off;
set(gcf,'Position',[100 100 800 300]);
 
figure;
plot(dates,h_hat);  xlim([dates(1) dates(end)]);  box off;
set(gcf,'Position',[100 100 800 300]);
 
j = randi(n); % or whatever other variable you fancy
figure;
plot(dates, [Y(:, j), F_hat'*Lambda_hat(j)]); xlim([dates(1) dates(end)]);  box off;
set(gcf,'Position',[100 100 800 300]);
 
toc;

