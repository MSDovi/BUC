function sv_block = SVBlock(dep_init, opt)
%%Function Description
% This function acts as a 'constructor' for an 'object' for a stochastic
% volatility block given by
%                  y_t = exp(0.5*h_t)e_t, e_t ~ N[0, 1]
%                  h_t = h_{t-1} + u_t, u_t ~ N[0, sigh2]    
%                  sigh2 ~ IG(nu, S)
%                  h_0 ~ N[a_0, b_0]
% dep refers to y_t, svsv (stochastic volatility state
% variable) refers to h_t, ic (initial condition) refers to h_0, var is
% sigh2.
% 
% If the option ncp is set to true, then a Normal prior for the standard
% deviation is considered, i.e.:
%
%                  sig ~ N[0, V_s]
%
% Inputs:
%   dep_init:     T x 1 vector containing an initialisation for the
%                 observed variable y of the stochastic volatility block
%   opt:          struct containing
%           nu:   nu parameter of IG prior for sigh2
%            S:   S parameter of IG prior for sigh2
%           Vs:   Prior variance for sigh in case of non-centred 
%                 parameterisation
%          ic0:   Prior mean for initial condition
%          icV:   Prior variance for initial condition
%          ncp:   Boolean indicating whether the stochastic volatility
%                 block should be estimated in noncentred 
%                 parameterisation
%
% Output:
%   sv_block:     struct containing
%        dep:     T x 1 vector of observed data (y_t in model above)
%        svsv:    T x 1 vector of stochastic volatiliy state variable 
%        ic:      Scalar value of initial value for stochastic volatility 
%                 state variable (svsv) 
%        var:     Scalar value of variance of volatility transition 
%                 equation drawn at the previous iteration (sigh2 in model 
%                 above)
%        T:       Scalar value saving dimension of the original dependent 
%                 variable.
%        opt:     struct passed in at construction
%        aux:     struct containing sparse differenced matrix HH

sv_block.opt = opt;
sv_block.dep = dep_init;
sv_block.T = length(sv_block.dep);
sv_block.ic = sv_block.opt.ic0;
if sv_block.opt.ncp
    sv_block.var = 0;
else
    sv_block.var = 1/(sv_block.opt.S * (sv_block.opt.nu - 1));
end

sv_block.svsv = log(sv_block.dep.^2 + .0001);
H = make_Hmat(sv_block.T, 1);
sv_block.aux.HH = H'*H;

if sv_block.opt.ncp
    sv_block.svsv_ncp = sv_block.svsv;
    sv_block.sd = (sv_block.var)^0.5;
    sv_block.sd_hat = 0;
    sv_block.Dsd = 0;

end