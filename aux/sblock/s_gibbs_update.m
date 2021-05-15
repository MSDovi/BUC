function s_block = s_gibbs_update(s_block)
%%Function Description
%
% This function samples s_t in the model described in SBlock.m
% 
% Input: 
%   s_block: SBlock 'object', see SBlock.m for a description.
% Output: 
%   s_block: SBlock 'object' with updated s, s0, alpha

%% Setup
y = s_block.dep;
T = s_block.T;
iSigma = inv(s_block.obsVar);
iOmega = inv(s_block.transVar);
H = s_block.H_coeff;
c0 = H\(s_block.alpha + s_block.exvar);
X = s_block.X;
iV_0 = inv(s_block.opt.icV);
ic_0 = s_block.opt.ic0;
p = s_block.p;

%% sample s
HiOgH = H'*iOmega*H;
Ks = HiOgH + iSigma;
s_hat = Ks\(HiOgH*c0 + iSigma*y);
s = s_hat + chol(Ks,'lower')'\randn(T,1);

%% sample s0
s_tilde = s - H\s_block.exvar;
Ks0 = X'*HiOgH*X + iV_0;
shat = Ks0\(X'*HiOgH*s_tilde + iV_0*ic_0);
s0 = shat + chol(Ks0,'lower')'\randn(p,1);

%% Update state block
s_block.s = s;
s_block.ic = s0;
s_block.alpha = s_block.X_tilde * s_block.ic;
s_block.res = s_block.H_coeff*s_block.s - s_block.alpha - s_block.exvar;

end