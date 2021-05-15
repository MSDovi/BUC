function df_block = DFBlock(dep_init, exvar_init, loadings_init, coeff_init, transVar_init, opt)
%%Function Description
% This function acts as a 'constructor' for an 'object' for a factor
% variable block given by the dual system:
%      y_t = \lambda f_t + e_t, e_t ~ N[0, Sigma_t], e ~ N[0, SIGMA]
%   f_t = mu_t + \Phi_1f_{t-1} + ... + \Phi_ps_{t-p} + u_t, u ~ N[0, Omega]
%                                                           u ~ N[0, OMEGA]
% y_t is a n x 1 vector of variables, \lambda is an n x q matrix of
% coefficients, f_t is a q x 1 vector, e_t is n x 1.^
%
% mu is some T x q vector of an exogenous variable (e.g., another state
% variable sampled previously).
%
% This object only samples f_t and \gamma = [f_0', ..., f_{-p+1}']'
%
% For generality, 'dep' refers to y_t, 'f' refers to f_t, ic 
% (initial conditions) refers to [f_0', f_{-1}', ..., f_{-p +1}']',
% 'obsVar' refers to Sigma, and 'transVar' refers to Omega.
% 
% Inputs:
%   dep_init:      T x n matrix containing an initialisation for the
%                  observed variable y 
%                  OR (class smart enough to realise which one you mean)
%                  Tn x 1 vector in the order [y_1', y_2', ..., y_T']'
%   exvar_init:    T x q vector containing an initialisation for mu
%                  OR (class smart enough to realise which one you mean)
%                  Tq x 1 vector in the order [mu_1', mu_2', ..., mu_T']'
%   loadings_init: n x q matrix of factor loadings
%   coeff_init:    q x q x p array containing an initialisation for the
%                  coefficients phi
%   obsVar_init:   T x T matrix for initalisation of observation equation
%                  OR (class smart enough to realise which one you mean)
%                  Tn x Tn matrix (in which case it should be sparse)
%   transVar_init: T x T matrix for initalisation of transition equation
%                  OR (class smart enough to realise which one you mean)
%                  Tq x Tq matrix (in which case it should be sparse)
%   opt:           struct containing
%                  ic0: pq x 1 vector containing prior means for initial
%                         conditions
%                       in the order [f_0', f_{-1}', ..., f_{-p + 1}]'
%                  icV: pq x pq covariance matrix for initial conditions
%                       in the order matching ic0
%                  ml:  (2*n-q-1)*q/2 * 1 vector of free factor loadings
%                  Vl:  (2*n-q-1)*q/2 x (2*n-q-1)*q/2 corresponding 
%                       covariance matrixmatrix of
%                  nu:  nu parameter for IG prior
%                  S:   S parameter for IG prior
%                  T:   T scalar (needed to help code figure out what
%                       notation you are passing in)
%
% Output:
%   df_block:        struct containing
%        dep:        Tn x 1 vector of observed data (y in model above)
%        Dep:        T x n matrix of reshaped dep
%        f:          Tq x 1 vector of factor variable
%        F:          T x q matrix of factor variable (exactly the same as f,
%                    just reshaped; useful for passing the output of this 
%                    block to another part in the Gibbs sampler)
%        exvar:      Tq x 1 vector for mu in transition equation
%        Exvar:      T x q matrix of reshaped exvar
%        loadings:   Tn x Tq matrix of factor loadings
%        loadings_s: n x T matrix of factor loadings
%        coeff:      q x q x p array of autoregressive coefficients
%        res:        Tq x 1 vector of residuals of the transition equation
%        Res:        T x q matrix of reshaped `res' residuals
%        Res_obs:    T x n matrix of residuals in observation equation
%        obsVar:     Tn x Tn covariance matrix for observation equation
%        transVar:   Tq x Tq covariance matrix for transition equation
%        H_coeff:    Tq x Tq matrix for representing autoregressive process
%        ic:         pq x 1 vector of initial conditions
%        F0:         p x q vector of initial conditions
%        alpha:      Tq x 1 vector containing initial conditions in matrix
%                    notation
%        X_tilde:    Tq x pq auxiliary matrix for sampling ic
%        X:          Tq x pq auxiliary matrix for sampling ic
%        p:          Scalar indicating number of autoregressive lags
%        n:          Scalar indicating number of elements in y_t
%        q:          Scalar indicating number of factors
%        T:          Scalar value indicating number of dimensions
%        opt:        struct passed in at construction

df_block.opt = opt;
df_block.coeff = coeff_init;
df_block.T = opt.T;
[junk, df_block.n] = size(dep_init);
[junk, df_block.q, df_block.p] = size(df_block.coeff);

obs_var_init = 1/(df_block.opt.S * (df_block.opt.nu - 1));
obsVar_init = sparse(1:df_block.n,1:df_block.n,(obs_var_init));

if size(dep_init, 1) == df_block.T
    df_block.dep = reshape(dep_init',df_block.T*df_block.n, 1);
else
    df_block.dep = dep_init;
end

if size(exvar_init, 1) == df_block.T
    df_block.exvar = reshape(exvar_init',df_block.T*df_block.q, 1);
else
    df_block.exvar = exvar_init;
end

if max(size(obsVar_init)) == df_block.n
    df_block.obsVar = kron(speye(df_block.T), obsVar_init);
elseif max(size(obsVar_init)) == (df_block.n * df_block.T)
    df_block.obsVar = obsVar_init;
end

if max(size(transVar_init)) == df_block.q
    df_block.transVar = kron(speye(df_block.T), transVar_init);
elseif max(size(transVar_init)) == (df_block.q * df_block.T)
    df_block.transVar = transVar_init;
end

df_block.H_coeff = df_make_Hmat(df_block.T, df_block.coeff);
df_block.ic = df_block.opt.ic0;
df_block.X_tilde = df_make_Xtilde(df_block.T, df_block.coeff);
df_block.alpha = df_block.X_tilde * df_block.ic;
df_block.X = df_block.H_coeff \ df_block.X_tilde;
df_block.f = df_block.exvar + df_block.X * df_block.ic;
df_block.res = df_block.H_coeff*df_block.f - df_block.alpha - df_block.exvar;
df_block.Res = reshape(df_block.res, df_block.q, df_block.T)';
df_block.loadings = kron(speye(df_block.T), loadings_init);
df_block.loadings_s = loadings_init;
df_block.F = reshape(df_block.f, df_block.q, df_block.T)';
df_block.F0 = reshape(df_block.ic, df_block.p, df_block.q);
df_block.Res_obs = reshape(df_block.dep - df_block.loadings*df_block.f, df_block.n, df_block.T)';

df_block.Dep = reshape(df_block.dep, df_block.n, df_block.T)';
df_block.Exvar = reshape(df_block.exvar, df_block.q, df_block.T)';


end