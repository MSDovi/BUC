function reg_block = RegBlock(dep_init, indie_init, opt)
%%Function Description
% Constructor for an object to draw from a standard Bayesian regression:
%               y = Xbeta + e, e ~ N[0, Sigma]
%                           beta ~ N[0, pVar]
%                          sigma ~ IG[nu, s]
%
% Inputs
%        dep_init:   T x 1 initialisation vector for y
%        indie_init: T x k initialisation matrix for X
%        opt:        struct containing
%                    pmean: k x 1 vector of prior means for beta
%                    pVar:  k x k matrix of prior covariance
%                    nu:    nu parameter for IG prior
%                    S:     S parameter for IG prior
%                     
% Output:
%   reg_block:       struct containing
%        dep:        T x 1 dependent vairable
%        indie:      T x k independent vairable
%        coeff:      k x 1 coefficient draw
%        res:        T x 1 vector of residuals
%        Var:        T x T variance draw
%        k:          Scalar indicating number of independent variables
%        T:          Scalar value indicating number of observations
%        opt:        struct passed in at construction
 
reg_block.opt = opt;
reg_block.dep = dep_init;
reg_block.T = length(reg_block.dep);
reg_block.k = max(size(opt.pmean));
reg_block.indie = indie_init;
reg_block.coeff = reg_block.opt.pmean;
%var_init = 1/(reg_block.opt.S * (reg_block.opt.nu - 1));
var_init = 1 / gamrnd(reg_block.opt.nu, reg_block.opt.S);
reg_block.Var = sparse(1:reg_block.T,1:reg_block.T,(var_init));
reg_block.res = reg_block.dep - reg_block.indie * reg_block.coeff;

end