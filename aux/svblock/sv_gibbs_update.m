function sv_block = sv_gibbs_update(sv_block)
%%Function Description
%
% This function samples h_t in the basic stochastic volatility model:
%
%                  y_t = exp(0.5*h_t)e_t, e_t ~ N[0, 1]
%                  h_t = h_{t-1} + u_t, u_t ~ N[0, sigh2]                       
%
% using a Gaussian mixture approximation as in Bayesian Econometric Methods
% (2nd edition) p. 391.
% 
% Input: 
%   sv_block: SVBlock 'object', see SVBlock.m for a description.
% Output: 
%   sv_block: SVBlock 'object' with updated svsv, ic, var

%% Setup
y = sv_block.dep;
h = sv_block.svsv;
h0 = sv_block.ic;
sigh2 = sv_block.var;
nu_h = sv_block.opt.nu;
S_h = sv_block.opt.S;
a0 = sv_block.opt.ic0;
b0 = sv_block.opt.icV;
HH = sv_block.aux.HH;
ystar = log(y.^2 + .0001);
T = sv_block.T;
if sv_block.opt.ncp
    Vs = sv_block.opt.Vs;
    sigh = sv_block.sd;
end

pj = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mj = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819]...
    - 1.2704;  
sigj2 = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
sigj = sqrt(sigj2);

%% Sample S from a 7-point distrete distribution
temprand = rand(T,1);
q = repmat(pj,T,1).*normpdf(repmat(ystar,1,7),...
    repmat(h,1,7)+repmat(mj,T,1),repmat(sigj,T,1));
q = q./repmat(sum(q,2),1,7);
S = 7 - sum(repmat(temprand,1,7)<cumsum(q,2),2) + 1;

%% Sample h, sigh2, h0 based on whether non-central parameterisation is used
if ~sv_block.opt.ncp

    %Sample h
    d_s = mj(S)'; 
    iSig_s = sparse(1:T,1:T,1./sigj2(S));
    Kh = HH/sigh2 + iSig_s;
    h_hat = Kh\(h0/sigh2*HH*ones(T,1) + iSig_s*(ystar-d_s));
    h = h_hat + chol(Kh,'lower')'\randn(T,1);

    %Sample sigh2
    sigh2 = 1/gamrnd(nu_h + T/2, 1/(S_h + (h-h0)'*HH*(h-h0)/2));
    
    % Sample h0
    Kh0 = 1/b0 + 1/sigh2;
    h0_hat = Kh0\(a0/b0 + h(1)/sigh2);
    h0 = h0_hat + sqrt(Kh0)'\randn;
    
else
    
    % Sample tilde h 
    d_s = mj(S)'; 
    iSig_s = sparse(1:T,1:T,1./sigj2(S));
    Kh = HH + sigh^2 * iSig_s;
    h_hat = Kh\(sigh*iSig_s*(ystar - d_s - h0));
    h_tilde = h_hat + chol(Kh,'lower')'\randn(T,1);
  
    % Sample sigh and h0
    Xbeta = [ones(T,1) h_tilde];
    iVbeta = diag([1/b0 1/Vs]);    
    Kbeta = iVbeta + Xbeta'*iSig_s*Xbeta;
    beta_hat = Kbeta\(iVbeta*[a0;0] + Xbeta'*iSig_s*(ystar-d_s));
    beta = beta_hat + chol(Kbeta,'lower')'\randn(2,1);
    h0 = beta(1); sigh = beta(2);
    
    % Randomly permute the signs h_tilde and sigh
    U = -1 + 2*(rand>0.5);
    h_tilde = U*h_tilde;
    sigh = U*sigh;
    sigh2 = sigh^2;

    % compute the mean and variance of the
    %    conditional density of omegah    
    Dbeta = Kbeta\speye(2);
    sigh_hat = beta_hat(2);
    Dsigh = Dbeta(2,2);
    
    h = h0 + sigh*h_tilde;    

end

%% Update volatility block
sv_block.svsv = h;
sv_block.var = sigh2;
sv_block.ic = h0;
if sv_block.opt.ncp
    sv_block.sd_hat = sigh_hat;
    sv_block.Dsd = Dsigh;
    sv_block.sd = sigh;
end

end