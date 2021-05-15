function ar_block = s_gibbs_update(ar_block)
%%Function Description
%
%This function samples coefficients in the model described in ARBlock.m
% 
%Input: 
%   s_block: ARBlock 'object', see ARBlock.m for a description.
%
%Output: 
%   s_block: ARBlock 'object' with new coefficients and variance draw

%% Setup
y = ar_block.dep - ar_block.exvar;
T = ar_block.T;
p = ar_block.p;
Xphi = ar_block.X_phi;
iSigma = inv(ar_block.Var);
phi0 = ar_block.opt.pmean;
iV = inv(ar_block.opt.pVar);

nu_sig2 = ar_block.opt.nu;
S_sig2 = ar_block.opt.S;

%% Sample coefficient vector
Kphi = iV + Xphi'*iSigma*Xphi;
phi_hat = Kphi\(iV*phi0 + Xphi'*iSigma*y);
phic = phi_hat + chol(Kphi, 'lower')\randn(p, 1);
if is_stat(phic)
    ar_block.countr = ar_block.countr + 1;
    ar_block.coeff = phic;    
end

%% Sample variance
var_draw = 1/gamrnd(nu_sig2 + T/2,1/(S_sig2 ...
        + (y-Xphi*ar_block.coeff)'*(y-Xphi*ar_block.coeff)/2));
ar_block.Var = sparse(1:T,1:T,(var_draw));

end