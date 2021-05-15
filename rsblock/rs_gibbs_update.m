function rs_block = rs_gibbs_update(rs_block)
%%Function Description
%
% This function samples the initial conditions for a 'residual' state
% variable. See RSBlock.m for a description of such state variables.
% 
% Input: 
%   rs_block: RSBlock 'object', see RSBlock.m for a description.
% Output: 
%   rs_block: RSBlock 'object' with updated ic and alpha.

%% Setup
iOmega = inv(rs_block.Var);
H = rs_block.H_coeff;
X = rs_block.X;
iV_0 = inv(rs_block.opt.icV);
ic_0 = rs_block.opt.ic0;
p = rs_block.p;
HiOgH = H'*iOmega*H;
s = rs_block.s;

%% sample s0
s_tilde = s - H\rs_block.exvar;
Ks0 = X'*HiOgH*X + iV_0;
shat = Ks0\(X'*HiOgH*s_tilde + iV_0*ic_0);
s0 = shat + chol(Ks0,'lower')'\randn(p,1);

%% Update state block
rs_block.ic = s0;
rs_block.alpha = rs_block.X_tilde * rs_block.ic;
rs_block.res = rs_block.H_coeff*rs_block.s - rs_block.alpha - rs_block.exvar;

end