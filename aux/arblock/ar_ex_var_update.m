function ar_block = ar_ex_var_update(ar_block, new_exvar, new_Var)
%%Function Description
%
% This function updates Sigma for objects of ARBlock class outside of the
% standard method included in reg_gibbs_update. This update is usually
% called when you want to pass in a covariance matrix computed through
% stochastic volatility.
%
% Inputs:
%   ar_block:    ARBlock 'object', see ARBlock.m for a description.
%   new_exvar:   T x 1 new exvar
%   new_Var:     T x T new variance object
%
% Output: 
%   ar_block:    ARBlock 'object' with updated covariance matrix

%% Update dependent variable

ar_block.exvar = new_exvar;
ar_block.Var = new_Var;

end