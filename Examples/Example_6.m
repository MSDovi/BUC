clear; clc;
tic;
%% Housekeeping
addpath('../aux')
addpath('../aux/svblock')
addpath('../aux/sblock')
%addpath('../aux/rsblock')
%addpath('../aux/arblock')
addpath('../aux/regblock')
intr_all = readmatrix('../data/interest.csv');
end_idx = 130;
intr_all = intr_all(1:end_idx, :);
intr_all(:, 1) = [];
pie_all = readmatrix('../data/inflation.csv');
pie_all = pie_all(1:end_idx, :);
pie_all(:, 1) = [];

countries = ["US", "Australia", "Belgium", "Canada", "France", "Germany", "Italy", "New_Zealand", "Norway", "Spain", "Switzerland", "UK"];

countries = ["France"]

nsim = 2000;
burnin = 500;
for i = 1:length(countries)
    countries(i)
    intr = intr_all(:, i);
    pi_e = pie_all(:, i);
    T = size(intr, 1);
    
    store_tau = zeros(nsim, T);
    store_tau_N = zeros(nsim, T);
    store_mu = zeros(nsim, 1);
    store_h = zeros(nsim, T);
    store_h_N = zeros(nsim, T);
    store_i_star = zeros(nsim, T);
    rej = zeros(nsim + burnin, T);

    %% Initialise the transition equation variance
    nu_sig2 = 0.1; S_sig2 = 1;
    sigma2_u = 1; 

    %% Initialise augmented data
    i_tilde_star = -abs(randn(T, 1));

    %% Initialise RegBlock
    opt_mu.pmean = 0; opt_mu.pVar = 0.1; % reasonable `economic' prior
    %opt_mu.pmean = min(intr); opt_mu.pVar = 0.000000001; % prior based on the fact that if they could have gone lower they would have
    opt_mu.nu = 0.1; opt_mu.S = T;
    reg_block = RegBlock(pi_e, ones(T, 1), opt_mu);

    %% Initialise SBlock
    opt_tau.ic0 = [mean(intr - pi_e)]; opt_tau.icV = speye(1);
    s_block_tau = SBlock(intr -  pi_e, zeros(T, 1), 1, ...
                       reg_block.Var, sigma2_u, opt_tau);

    %% Initialise Naive SBlock
    nu_N = 0.1; S_N = 1;
    sigma2_epsN = 3;
    sigma2_uN = 3;
    s_block_tau_N = SBlock(intr - pi_e, zeros(T, 1), 1, ...
                       reg_block.Var, sigma2_uN, opt_tau);
    %% Initialise SVBlocks
    opt_sv.nu = 3; opt_sv.S = T; opt_sv.Vs = 0.1;
    opt_sv.ic0 = 1; opt_sv.icV = 100;
    opt_sv.ncp = true;
    sv_block_h = SVBlock(intr - s_block_tau.s, opt_sv);
    sv_block_h_N = SVBlock(intr - s_block_tau_N.s, opt_sv);
    sv_block_tau = SVBlock(s_block_tau.res, opt_sv);
    %% Gibbs Sampler
    
    for isim = 1:nsim + burnin
        %Sigma = sparse(1:T,1:T,exp(sv_block_h.svsv));
        Sigma = reg_block.Var;
        %Sigma_N = s_block_N.Var;
        Sigma_N = sparse(1:T,1:T,exp(sv_block_h_N.svsv));  
        Omega = sparse(1:T,1:T,(sigma2_uN));  
        %Omega = sparse(1:T,1:T,exp(sv_block_tau.svsv));  
        mu = reg_block.coeff;
        i_tilde = intr - mu;
        i_tilde_z = (i_tilde <= 0).*i_tilde_star + (i_tilde > 0).*i_tilde;

        % Sample \tau and \tau_0
        i_hat_tau = i_tilde_z - pi_e + mu;
        s_block_tau = s_ex_coeff_var_update(s_block_tau, zeros(T, 1), 1, Sigma, Omega);
        s_block_tau = s_dep_update(s_block_tau, i_hat_tau);
        s_block_tau = s_gibbs_update(s_block_tau);

        % Sample \sigma^2_u
        sigma2_u = 1/gamrnd(nu_sig2 + T/2,1/(S_sig2 ...
            + s_block_tau.res'*s_block_tau.res/2));

        % Sample \mu^2 and \sigma^_\varepsilon
        tau = s_block_tau.s;
        i_hat_mu = i_tilde_z - pi_e - tau;
        reg_block = reg_dep_indie_update(reg_block, -1*i_hat_mu, ones(T, 1));
        reg_block = reg_var_update(reg_block, Sigma);
        reg_block = reg_gibbs_update(reg_block);    

        % Sample i_tilde_star
        for t = 1:T
           cand = -mu + pi_e(t) + tau(t) + sqrt(Sigma(t,t))*randn;
           if cand <= 0
            i_tilde_star(t) = cand;
            rej(isim, t) = rej(isim, t) + 1;
           end
        end

        % Naive trend
        s_block_tau_N = s_ex_coeff_var_update(s_block_tau_N, zeros(T, 1), 1, Sigma_N, sparse(1:T, 1:T, sigma2_uN));
        s_block_tau_N = s_gibbs_update(s_block_tau_N);
        sigma2_uN = 1/gamrnd(nu_N + T/2,1/(S_N ...
            + s_block_tau_N.res'*s_block_tau_N.res/2));   
        sigma2_epsN = 1/gamrnd(nu_sig2 + T/2,1/(S_sig2 ...
            + (s_block_tau_N.dep - s_block_tau_N.s)'*(s_block_tau_N.dep - s_block_tau_N.s))); 
        
        % Stochastic Volatility
        sv_block_h = sv_dep_update(sv_block_h, i_tilde_z + mu - pi_e - tau);
        sv_block_h = sv_gibbs_update(sv_block_h);    
        sv_block_h_N = sv_dep_update(sv_block_h_N, s_block_tau_N.res);
        sv_block_h_N = sv_gibbs_update(sv_block_h_N);      
        sv_block_tau = sv_dep_update(sv_block_tau, s_block_tau.res);
        sv_block_tau = sv_gibbs_update(sv_block_tau);      
        
        if isim > burnin
                isave = isim - burnin;
                store_tau(isave, :) = tau';
                store_tau_N(isave, :) = s_block_tau_N.s';
                store_mu(isave) = mu;
                store_i_star(isave, :) = (i_tilde_star + mu)';
                store_h_N(isave, :) = sv_block_h_N.svsv';
                store_h(isave, :) = sv_block_h.svsv';
        end

         if (mod(isim,2000) == 0)
            disp([num2str(isim) ' loops... '])
        end

    end

    tau_hat = mean(store_tau)';
    tau_hat_N = mean(store_tau_N)';
    mu_hat = mean(store_mu);
    istar_hat = mean(store_i_star)';
    h_hat = mean(store_h)';
    h_N_hat = mean(store_h_N)';

    %dates = datetime(1954, 7, 1) + calmonths(1:T);
    dates = linspace(1980.25, 2019.5, T);
    figure; hold on
    a1 = plot(dates, intr-pi_e); M1 = "Real Interest Rate";
    a2 = plot(dates, tau_hat); M2 = "Trend";
    a3 = plot(dates, tau_hat_N); M3 = "Naive Trend";
    legend([a1, a2, a3], [M1, M2, M3]);
%    saveas(gcf, strcat('ZLBFigs/', countries(i), 'Trends.png'))

    figure; hold on
    a2 = plot(dates, tau_hat - tau_hat_N); M2 = "Trend - Naive trend";
    legend([a2], [M2]);
%    saveas(gcf, strcat('ZLBFigs/', countries(i), 'TrendDiff.png'))

%     figure; hold on
%     a1 = plot(dates, h_hat); M1 = "Stochastic Volatility";
%     a2 = plot(dates, h_N_hat); M2 = "Stochastic Volatility Naive";
%     legend([a1, a2], [M1, M2]);
%     saveas(gcf, strcat('ZLBFigs/', countries(i), 'SV.png'))
    
    figure; hold on
    a2 = plot(dates, tau_hat - tau_hat_N); M2 = "Trend - Naive trend";
    legend([a2], [M2]);
%    saveas(gcf, strcat('ZLBFigs/', countries(i), 'TrendDiff.png'))
    
%     figure; hold on
%     a2 = plot(dates, tau_hat); M2 = "Trend";
%     legend([a2], [M2]);
%     saveas(gcf, strcat('ZLBFigs/', countries(i), 'TrendSOLO.png'))
    
    mu_hat
end
