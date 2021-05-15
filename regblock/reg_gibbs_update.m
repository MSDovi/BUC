function reg_block = reg_gibbs_update(reg_block)
%%Function Description
%
% This function samples coefficients in the model described in RegBlock.m
% 
% Input: 
%   s_block: RegBlock 'object', see RegBlock.m for a description.
% Output: 
%   s_block: RegBlock 'object' with new coefficients and variance draw

%% Setup
y = reg_block.dep;
T = reg_block.T;
k = reg_block.k;
X= reg_block.indie;
Sigma = reg_block.Var;
beta0 = reg_block.opt.pmean;
iV = inv(reg_block.opt.pVar);
nu_sig2 = reg_block.opt.nu;
S_sig2 = reg_block.opt.S;

%% Sample coefficient vector
Kbeta = iV + X'*(Sigma\X);
beta_hat = Kbeta\(iV*beta0 + X'*(Sigma\y));
beta = beta_hat + chol(Kbeta, 'lower')\randn(k, 1);

%% Sample variance
var_draw = 1/gamrnd(nu_sig2 + T/2,1/(S_sig2 ...
        + (y-X*beta)'*(y-X*beta)/2));

%% Save draws
reg_block.coeff = beta;
reg_block.Var = sparse(1:T,1:T,(var_draw));
reg_block.res = reg_block.dep - reg_block.indie * reg_block.coeff;

end