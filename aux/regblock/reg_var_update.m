function reg_block = reg_var_update(reg_block, new_Var)
%%Function Description
%
% This function updates Sigma for objects of RegBlock class outside of the
% standard method included in reg_gibbs_update. This update is usually
% called when you want to pass in a covariance matrix computed through
% stochastic volatility.
% Inputs:
%   reg_block: RegBlock 'object', see RegBlock.m for a description.
%   new_Sigma:  T x T new covariance matrix
% Output: 
%   s_block: RegBlock 'object' with updated covariance matrix

%% Update dependent variable
reg_block.Var = new_Var;
end