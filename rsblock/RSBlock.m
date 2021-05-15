function rs_block = RSBlock(exvar_init, coeff_init, Var_init, opt)
%%Function Description
% This function acts as a 'constructor' for an 'object' for a  'residual' 
% state variable, i.e. a state variable that is perfectly determined
% conditional on the other parts of the UC model. This object is
% essentially only used to sample the initial conditions of this 'residual'
% state vairable
%
% The block is given by the single unobserved transition equation of the
% residual state variable, s:
%   s_t = mu_t + \phi_1s_{t-1} + ... + \phi_ps_{t-p} + u_t, u ~ N[0, Omega]
% 's' refers to s_t, ic (initial conditions) refers to [s_0, s_{-1}, ..., 
% s_{-p +1}], and 'Var' refers to Omega.
% Inputs:
%   coeff_init:    p x 1 vector containing an initialisation for the
%                  coefficients phi
%   exvar_init:    T x 1 vector containing an initialisation for mu
%   Var_init:      T x T matrix for initalisation of transition equation
%   opt:           struct containing
%                  ic0: p x 1 vector containing prior means for initial
%                         conditions
%                  icV: p x p covariance matrix for initial conditions
% Output:
%   s_block:       struct containing
%        s:        T x 1 vector of state variable
%        exvar:    T x 1 vector for mu 
%        coeff:    p x 1 vector of autoregressive coefficients
%        res:      T x 1 vector of residuals of the transition equation
%        Omega:    T x T covariance matrix for transition equation
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

rs_block.opt = opt;
rs_block.T = length(exvar_init);
rs_block.coeff = coeff_init;
rs_block.H_coeff = make_Hmat(rs_block.T, rs_block.coeff);
rs_block.p = length(rs_block.coeff);
rs_block.ic = rs_block.opt.ic0;
rs_block.X_tilde = rs_make_Xtilde(rs_block.T, rs_block.coeff);
rs_block.alpha = rs_block.X_tilde * rs_block.ic;
rs_block.X = rs_block.H_coeff \ rs_block.X_tilde;
rs_block.Var = Var_init;
rs_block.exvar = exvar_init;
rs_block.s = rs_block.exvar + rs_block.X * rs_block.ic;
rs_block.res = rs_block.H_coeff*rs_block.s - rs_block.alpha - rs_block.exvar;

end