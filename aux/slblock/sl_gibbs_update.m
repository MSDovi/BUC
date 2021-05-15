function sl_block = sl_gibbs_update(sl_block)
%%Function Description
%
% This function samples s_t in the model described in SLBlock.m
% 

%% Setup
y = sl_block.dep - sl_block.xi;
T = sl_block.T;
iSigma = inv(sl_block.obsVar);
iOmega = inv(sl_block.transVar);
H = sl_block.H_coeff;
G = sl_block.G_coeff;
c0 = H\(sl_block.alpha + sl_block.exvar);
X = sl_block.X;
iV_0 = inv(sl_block.opt.icV);
ic_0 = sl_block.opt.ic0;
q = length(ic_0);

%% sample s
HiOgH = H'*iOmega*H;
Ks = HiOgH + G'*iSigma*G;
s_hat = Ks\(HiOgH*c0 + G'*iSigma*y);
s = s_hat + chol(Ks,'lower')'\randn(T,1);

%% sample s0
s_tilde = s - H\sl_block.exvar;
Ks0 = X'*HiOgH*X + iV_0;
shat = Ks0\(X'*HiOgH*s_tilde + iV_0*ic_0);
s0 = shat + chol(Ks0,'lower')'\randn(q,1);

%% Update state block
sl_block.s = s;
sl_block.ic = s0;
ic_lags_obs = sl_block.ic(1:sl_block.p_1)
ic_leads = sl_block.ic(sl_block.q_1 + 1:sl_block.q_1 + sl_block.p_2)
sl_block.ic_obs = [ic_lags_obs'; ic_leads'];
sl_block.alpha = sl_block.X_tilde * sl_block.ic;
sl_block.xi = sl_block.X_tilde_obs * sl_block.ic_obs;
sl_block.res = sl_block.H_coeff*sl_block.s - sl_block.alpha - sl_block.exvar;

end