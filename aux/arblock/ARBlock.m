function ar_block = ARBlock(dep_init, exvar_init, ic_init, opt)
%%Function Description
% Constructor for an object to takes draws from a stationary Bayesian 
% autoregressive process:
%     y_t = mu_t + phi_1y_{t-1} + ... + phi_py_{t-p} + e_t, e ~ N[0, Sigma]
%
% mu is some T x 1 vector of an exogenous variable (e.g., a state
% variable sampled previously). 
%
% Inputs
%        dep_init:   T x 1 initialisation vector
%        exvar_init: T x 1 vector containing an initialisation for mu
%        ic_init:    p x 1 initialisation vector for initial conditions
%        opt:        struct containing
%                    pmean: p x 1 vector of prior means
%                    pVar:  p x p matrix of prior covariance
%                    nu:    nu parameter for IG prior
%                    S:     S parameter for IG prior
%                     
% Output:
%    ar_block:       struct containing
%        dep:        T x 1 dependent vairable
%        coeff:      p x 1 coefficient draw
%        Var:        T x T covariance matrix
%        X_phi:      T x p matrix for sampling the phi coefficients 
%                    (the coefficients are sampled in another block)
%        p:          Scalar indicating number of lags
%        T:          Scalar value indicating number of observations
%        opt:        struct passed in at construction

ar_block.opt = opt;

ar_block.dep = dep_init;
ar_block.T = length(ar_block.dep);
ar_block.p = max(size(opt.pmean));
ar_block.coeff = ar_block.opt.pmean;
ar_block.ic = ic_init;
ar_block.X_phi = ar_make_Xphi(ar_block.dep, ar_block.ic);
ar_block.countr = 0;
ar_block.exvar = exvar_init;
var_init = 1/(ar_block.opt.S * (ar_block.opt.nu - 1));
ar_block.Var = sparse(1:ar_block.T,1:ar_block.T,(var_init));

end