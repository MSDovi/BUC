function s_block = s_ex_coeff_var_update(s_block, new_exvar, new_coeff, new_obsVar, new_transVar)
%%Function Description
% This function updates the fields mu, phi, Sigma, and Omega for an object of
% the class SBlock (see SBlock.m for details)

s_block.exvar = new_exvar;
s_block.coeff = new_coeff;
s_block.obsVar = new_obsVar;
s_block.transVar = new_transVar;

s_block.H_coeff = make_Hmat(s_block.T, s_block.coeff);
s_block.X_tilde = s_make_Xtilde(s_block.T, s_block.coeff);
s_block.X = s_block.H_coeff \ s_block.X_tilde;

end