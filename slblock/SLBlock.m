function sl_block = SLBlock(dep_init, exvar_init, coeff_init, obsVar_init, transVar_init, opt)
%%Function Description
% This 'class' is a generalisation of the SBlock class. In particular, it
% considers the model
%		y_t &= cs_t + a_1s_{t-1} + a_2s_{t-2} + \dots + a_{p_1}s_{t-p_1} + b_1s_{t+1} + b_2s_{t+2} + \dots + b_{p_2}s_{t+p_2} + \varepsilon_t, \varepsilon \sim \mathcal{N}[0, \Sigma]\\
%		s_t &= \mu_t + \phi_1s_{t-1} + \phi_2s_{t-2} + \dots + \phi_{q_1}s_{t-q_1} + \theta_1s_{t+1} + \theta_2s_{t+2} + \dots + \theta_{q_2}s_{t + q_2} + v_t, v_t \sim \mathcal{N}[0, \Omega]%
%
%
% (it may be more readable to look at this in the documentation folder). 
%
% mu is some T x 1 vector of an exogenous variable (e.g., another state
% variable sampled previously).
%
% For generality, 'dep' refers to y_t, 's' refers to s_t, svsv, ic 
% (initial conditions) refers to $[s_0, s_{-1}, \dots, s_{-q_1 + 1}, s_{T + 1}, s_{T + 2}, \dots, S_{T + q_2}]'$, 
% 'obsVar' refers to Sigma, and 'transVar' refers to Omega.
% Inputs:
%   dep_init:      T x 1 vector containing an initialisation for the
%                  observed variable y 
%   exvar_init:    T x 1 vector containing an initialisation for mu
%   coeff_init:    struct containing
%                  zero: T x 1 coefficient for c in equation above
%                  obs_lag: p_1 x 1 vector of coefficients for a
%                  obs_lead: p_2 x 1 vector of coefficients for b
%                  trans_lag: q_1 x 1 vector of coefficients for \phi
%                  trans_lead: q_2 x 1 vector of coefficients for \theta
%   obsVar_init:   T x T matrix for initalisation of observation equation
%   transVar_init: T x T matrix for initalisation of transition equation
%   opt:           struct containing
%                  ic0: (q_1 + q_2) vector containing prior means for initial
%                         conditions
%                  icV: (q_1 + q_2) covariance matrix for initial conditions
%
% Output:
%   s_block:       struct containing
%        dep:      T x 1 vector of observed data (y in model above)
%        s:        T x 1 vector of state variable
%        exvar:    T x 1 vector for mu in transition equation
%        coeff:    struct containing
%                  zero: scalar coefficient for c in equation above
%                  obs_lag: p_1 x 1 vector of coefficients for a
%                  obs_lead: p_2 x 1 vector of coefficients for b
%                  trans_lag: q_1 x 1 vector of coefficients for \phi
%                  trans_lead: q_2 x 1 vector of coefficients for \theta
%        res:      T x 1 vector of residuals of the transition equation
%        obsVar:   T x T covariance matrix for observation equation
%        transVar: T x T covariance matrix for transition equation
%        G_coeff:  T x T auxiliary matrix for observation equation
%        H_coeff:  T x T matrix for representing autoregressive process
%        ic:       (q_1 + q_2) x 1 vector of initial and final conditions
%        ic_obs:   (p_1 + p_2) x 1 vector of initial and final conditions
%                  for observation equation; these will be a subset of the
%                  ones in ic.
%        alpha:    T x 1 vector containing initial conditions in matrix
%                  notation
%        X_tilde:  T x p auxiliary matrix for sampling ic
%        X:        T x p auxiliary matrix for sampling ic
%        T:        Scalar value saving dimension of the original dependent 
%                  variable.
%        opt:      struct passed in at construction

sl_block.opt = opt;
sl_block.dep = dep_init;
sl_block.obsVar = obsVar_init;
sl_block.T = length(sl_block.dep);
sl_block.exvar = exvar_init;
sl_block.coeff = coeff_init;

sl_block.q_1 = length(sl_block.coeff.trans_lag);
sl_block.q_2 = length(sl_block.coeff.trans_lead);
sl_block.p_1 = length(sl_block.coeff.obs_lag);
sl_block.p_2 = length(sl_block.coeff.obs_lead);

sl_block.H_coeff = sl_make_Hmat(sl_block.T, sl_block.coeff.trans_lag, sl_block.coeff.trans_lead, 1);
sl_block.G_coeff = sl_make_Hmat(sl_block.T, -1*sl_block.coeff.obs_lag, -1*sl_block.coeff.obs_lead, sl_block.coeff.zero);
sl_block.ic = sl_block.opt.ic0;
ic_lags_obs = sl_block.ic(1:sl_block.p_1);
ic_leads = sl_block.ic(sl_block.q_1 + 1:sl_block.q_1 + sl_block.p_2);
sl_block.ic_obs = [ic_lags_obs; ic_leads];

sl_block.X_tilde = sl_make_Xtilde(sl_block.T, sl_block.coeff.trans_lag, sl_block.coeff.trans_lead);
sl_block.X_tilde_obs = sl_make_Xtilde(sl_block.T, sl_block.coeff.obs_lag, sl_block.coeff.obs_lead);
sl_block.alpha = sl_block.X_tilde * sl_block.ic;
sl_block.xi = sl_block.X_tilde_obs * sl_block.ic_obs;

sl_block.X = sl_block.H_coeff \ sl_block.X_tilde;
sl_block.s = sl_block.exvar + sl_block.X * sl_block.ic;
sl_block.transVar = transVar_init;
sl_block.res = sl_block.H_coeff*sl_block.s - sl_block.alpha - sl_block.exvar;

end