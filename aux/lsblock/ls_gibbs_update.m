 function ls_block = ls_gibbs_update(ls_block)
 
    sigtau = ls_block.sig_trend;
    sigalp = ls_block.sig_drift;
    sigy2 = ls_block.sig_obs;
    T = ls_block.T;
    Pgam = ls_block.Pgam;
    alp0 = ls_block.alp0;
    tau0 = ls_block.tau0;
    y = ls_block.dep;
    beta0 = ls_block.beta0;
    iVbeta = ls_block.iVbeta;
    H = ls_block.H;
    nu_sig0 = ls_block.nu_sig0;
    S_sig0 = ls_block.S;
    

    Xgam = [sigtau*speye(T) sigalp*speye(T)];   
    Kgam = Pgam + Xgam'*Xgam/sigy2;    
    gam_hat = Kgam\(1/sigy2*Xgam'*(y-tau0-(1:T)'*alp0));    
    gam = gam_hat + chol(Kgam,'lower')'\randn(2*T,1); 
        
    % sample beta
    X = [ones(T,1) (1:T)' gam(1:T) gam(T+1:end)];
    beta = ls_sample_coeffs(y,X,sigy2,beta0,iVbeta);
    tau0 = beta(1); alp0 = beta(2); sigtau = beta(3); sigalp = beta(4);
    
    % permutate the signs of gam and (sigtau, sigalp)        
    U = -1 + 2*(rand>0.5);
    gam = U*gam;
    sigalp = U*sigalp;
    sigtau = U*sigtau;
    
    % compute tau and alp
    tau_tilde = gam(1:T);
    A_tilde = gam(T+1:end);
    alp_tilde = H*A_tilde;
    alp = alp0 + sigalp*alp_tilde;
    tau = tau0 + sigtau*tau_tilde + (1:T)'*alp0 + sigalp*(H\alp_tilde);
    
    sigy2 = 1/gamrnd(nu_sig0 + T/2,1/(S_sig0 + (y-tau)'*(y-tau)/2)); 
    
    ls_block.trend = tau;
    ls_block.drift = alp;
    ls_block.sig_obs = sigy2;
    ls_block.alp0 = alp0;
    ls_block.tau0 = tau0;    
    ls_block.sig_trend = sigtau;
    ls_block.sig_drift = sigalp;

       
 end
