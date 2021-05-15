function s_block = SBlock(dep_init, exvar_init, coeff_init, obsVar_init, transVar_init, opt)
%%Function Description
% This function acts as a 'constructor' for an 'object' for a state variable
% block given by the dual system:
%      y_t = c_t + e_t, e ~ N[0, Sigma]
%   s_t = mu_t + \phi_1s_{t-1} + ... + \phi_ps_{t-p} + u_t, u ~ N[0, Omega]
% 
% mu is some T x 1 vector of an exogenous variable (e.g., another state
% variable sampled previously).
%
% For generality, 'dep' refers to y_t, 's' refers to s_t, svsv, ic 
% (initial conditions) refers to [s_0, s_{-1}, ..., s_{-p +1}], 'obsVar' 
% refers to Sigma, and 'transVar' refers to Omega.
% Inputs:
%   dep_init:      T x 1 vector containing an initialisation for the
%                  observed variable y 
%   exvar_init:    T x 1 vector containing an initialisation for mu
%   coeff_init:    p x 1 vector containing an initialisation for the
%                  coefficients phi
%   obsVar_init:   T x T matrix for initalisation of observation equation
%   transVar_init: T x T matrix for initalisation of transition equation
%   opt:           struct containing
%                  ic0: p x 1 vector containing prior means for initial
%                         conditions
%                  icV: p x p covariance matrix for initial conditions
%
% Output:
%   s_block:       struct containing
%        dep:      T x 1 vector of observed data (y in model above)
%        s:        T x 1 vector of state variable
%        exvar:    T x 1 vector for mu in transition equation
%        coeff:    p x 1 vector of autoregressive coefficients
%        res:      T x 1 vector of residuals of the transition equation
%        obsVar:   T x T covariance matrix for observation equation
%        transVar: T x T covariance matrix for transition equation
%        H_coeff:  T x T matrix for representing autoregressive process
%        ic:       p x 1 vector of initial conditions
%        alpha:    T x 1 vector containing initial conditions in matrix
%                  notation
%        X_tilde:  T x p auxiliary matrix for sampling ic
%        X:        T x p auxiliary matrix for sampling ic
%        p:        Scalar indicating number of autoregressive lags
%        T:        Scalar value saving dimension of the original dependent 
%                  variable.
%        opt:      struct passed in at construction

s_block.opt = opt;
s_block.dep = dep_init;
s_block.obsVar = obsVar_init;
s_block.T = length(s_block.dep);
s_block.exvar = exvar_init;
s_block.coeff = coeff_init;

s_block.H_coeff = make_Hmat(s_block.T, s_block.coeff);
s_block.p = length(s_block.coeff);

s_block.ic = s_block.opt.ic0;
s_block.X_tilde = s_make_Xtilde(s_block.T, s_block.coeff);
s_block.alpha = s_block.X_tilde * s_block.ic;
s_block.X = s_block.H_coeff \ s_block.X_tilde;
s_block.s = s_block.exvar + s_block.X * s_block.ic;
s_block.transVar = transVar_init;
s_block.res = s_block.H_coeff*s_block.s - s_block.alpha - s_block.exvar;

end