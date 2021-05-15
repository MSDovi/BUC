%% Housekeeping
clear; clc;
tic;
addpath('../aux')
addpath('../aux/svblock')
addpath('../aux/sblock')
addpath('../aux/rsblock')
addpath('../aux/arblock')
addpath('../aux/dfblock')

raw = xlsread('../data/cpi_core_sa_OECD.xlsx');
raw = raw(81:end, :); % remove NAs

data = 400 * log( raw( 2 : end , : ) ./ raw( 1 : end - 1, : ) ) ;


Y = data;
[T, n] = size(Y);
y = reshape(Y',T*n,1);




cpi = readmatrix('../data/China/CPIHolzSun.csv');

provinces = ["Beijing", "Anhui", "Chongqing", "Fujian", "Gansu", "Guangdong", ...
             "Guangxi",  "Guizhou", "Hainan", "Hebei",  "Heilongjiang", ...
             "Henan", "Hubei" "Hunan" "InnerMongolia", "Jiangsu", "Jiangxi", ...
             "Jilin", "Liaoning", "Ningxia", "Qinghai", "Shaanxi",  "Shandong", ...
             "Shanghai", "Shanxi", "Sichuan", "Tianjin", "Tibet", "Xinjiang", ...
             "Yunnan", "Zhejiang"];
raw = cpi(:, 2:end);
years = cpi(2:end, 1);
Y = 100 * log( raw( 2 : end , : ) ./ raw( 1 : end - 1, : ) ) ;
[T, n] = size(Y);
y = reshape(Y',T*n,1);













%data_mean = mean(data)';
%data_std = std(data)';
%idx = ~isnan(data_mean); % extract variables that don't have missing data
%Y = (data(:,idx) - repmat(data_mean(idx)',T,1))./repmat(data_std(idx)',T,1); % standardised the data
 
nsim = 1000;
burnin = 1000;
 
store_F = zeros(nsim,T);
store_Lambda = zeros(nsim,n);
store_c = zeros(nsim, T, n);
store_tau = zeros(nsim, T, n);
store_g = zeros(nsim,T);
store_h = zeros(nsim, T, n);
store_l = zeros(nsim, T, n);

%% Initialise DFBlock
nLambda = (2*n-2)/2;
opt.ic0 = 0; opt.icV = 0.5;
opt.ml = zeros(nLambda, 1); opt.Vl = 100 * speye(nLambda);
opt.nu = 5; opt.S = 1*(opt.nu-1);
opt.T = T;
phi_f = 0.6;
Lambda_init = zeros(n,1);
Lambda_init(1:n+1:end) = 1; % set diagonal elements to 1

df_block = DFBlock(Y, zeros(T, 1), Lambda_init, phi_f, sparse(diag(1)), opt);

%% Initialise ARBlocks
opt_phi.pmean = [0.6]; opt_phi.pVar = speye(1);
opt_phi.nu = 3; opt_phi.S = 1*(opt_phi.nu-1);

for i = 1:n
    ar_block_c(i) = ARBlock(0.01*Y(:, i), zeros(T, 1), zeros(1, 1), opt_phi);
end

ar_block_f = ARBlock(df_block.f, zeros(T, 1), zeros(1, 1), opt_phi);

%% Initialise State Variables Block
opt_c.ic0 = [0]'; opt_c.icV = 10*speye(1);
for i = 1:n
    s_block_c(i) = SBlock(0.01*Y(:, i), zeros(T, 1), opt_phi.pmean, sparse(diag(1)), ar_block_c(i).Var, opt_c);
end
               
%% Initialise the Residual State Variables Block
 opt_tau.icV = speye(1);
 for i = 1:n
     opt_tau.ic0 = [0];
     rs_block_tau(i) = RSBlock(zeros(T, 1), 1, sparse(diag(1)), opt_tau);
 end

 %% Initialise SVBlocks
opt_sv.nu = 3; opt_sv.S = 0.2; opt_sv.Vs = 0.1;
opt_sv.ic0 = 0; opt_sv.icV = 100;
opt_sv.ncp = true;
 
sv_block_g = SVBlock(df_block.f - df_block.exvar, opt_sv);
 
for i = 1:n
    sv_block_h(i) = SVBlock(rs_block_tau(i).res, opt_sv);
end

for i = 1:n
    sv_block_l(i) = SVBlock(s_block_c(i).res, opt_sv);
end

%% Initialise some things defined purely for readability
flat_Sigma = zeros(T*n, 1);
flat_Omega = zeros(T*n, 1);
Tau = Y;
Taum1 = Tau;
C = Y - Tau;

%% Gibbs Sampler
for isim = 1:nsim + burnin
    isim
    Xi = sparse(1:T,1:T,exp(sv_block_g.svsv)); % T x T
    for i = 1:n
       slct = i:n:((T-1)*n+i);
       flat_Sigma(slct) = exp(sv_block_h(i).svsv);
       flat_Omega(slct) = exp(sv_block_l(i).svsv);
    end
    Sigma = sparse(1:(T*n), 1:(T*n), flat_Sigma); % Tn x Tn
    Omega = sparse(1:(T*n), 1:(T*n), flat_Omega); % Tn x Tn

    phi_f = ar_block_f.coeff;
        
    % sample f and factor loadings
    df_block = df_ex_coeff_var_update(df_block, zeros(T, 1), phi_f, Sigma, Xi);
    df_block = df_dep_update(df_block, Y - Taum1 - C);
    df_block = df_gibbs_update(df_block);
    % Unpack some things for easier reading
    F = df_block.F; F_0 = df_block.F0;
    Lambda = df_block.loadings_s;
    
    % Sample c, c_0 (and back out Tau), and sample Tau0 (and back out
    % Taum1)
    for i = 1:n
        slct = i:n:((T-1)*n+i);
        Sigma_j = sparse(1:T, 1:T, flat_Sigma(slct));
        Omega_j = sparse(1:T, 1:T, flat_Omega(slct));
        s_block_c(i) = s_ex_coeff_var_update(s_block_c(i), zeros(T, 1), ar_block_c(i).coeff, Sigma_j, Omega_j);
        s_block_c(i) = s_dep_update(s_block_c(i), Y(:, i) - Lambda(i)*F - Taum1(:, i));
        s_block_c(i) = s_gibbs_update(s_block_c(i));
        Tau(:, i) = Y(:, i) - Lambda(i)*F - (s_block_c(i).s);
        rs_block_tau(i) = rs_ex_coeff_var_update(rs_block_tau(i), zeros(T, 1), 1, Sigma_j);
        rs_block_tau(i) = rs_s_update(rs_block_tau(i), Tau(:, i));
        rs_block_tau(i) = rs_gibbs_update(rs_block_tau(i));
        Taum1(:, i) = [rs_block_tau(i).ic; Tau(1:end-1, i)];
        
        ar_block_c(i) = ar_dep_ic_update(ar_block_c(i), s_block_c(i).s, s_block_c(i).ic);
        ar_block_c(i) = ar_ex_var_update(ar_block_c(i), zeros(T, 1), Omega_j);
        ar_block_c(i) = ar_gibbs_update(ar_block_c(i));
        
    end
    
    % sample Phis
    ar_block_f = ar_dep_ic_update(ar_block_f, df_block.f, df_block.ic);
    ar_block_f = ar_ex_var_update(ar_block_f, zeros(T, 1), Xi);
    ar_block_f = ar_gibbs_update(ar_block_f);
 
    % SV Block
    for i = 1:n
        sv_block_h(i) = sv_dep_update(sv_block_h(i), rs_block_tau(i).res);
        sv_block_h(i) = sv_gibbs_update(sv_block_h(i));
        
        sv_block_l(i) = sv_dep_update(sv_block_l(i), s_block_c(i).res);
        sv_block_l(i) = sv_gibbs_update(sv_block_l(i));     
    end
    sv_block_g = sv_dep_update(sv_block_g, df_block.Res);
    sv_block_g = sv_gibbs_update(sv_block_g);
    
    if isim > burnin
        isave = isim - burnin;
        store_F(isave,:,:) = F;
        store_Lambda(isave,:,:) = Lambda;        
        store_g(isave,:) = sv_block_g.svsv'; 
        for i = 1:n
            store_h(isave, :, i) = sv_block_h(i).svsv';
            store_l(isave, :, i) = sv_block_l(i).svsv';
            store_c(isave, :, i) = s_block_c(i).s;
            store_tau(isave, :, i) = rs_block_tau(i).s;
        end
    end
    
     if (mod(isim,2000) == 0)
        disp([num2str(isim) ' loops... '])
    end
    
end

tau_hat = squeeze(mean(store_tau));
c_hat = squeeze(mean(store_c));

figure;
plot([tau_hat(:, 1), Y(:, 1)]);

figure;
plot(c_hat(:, 1));

figure;
plot(Y(:, 1));




