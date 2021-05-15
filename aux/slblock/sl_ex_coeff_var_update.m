function sl_block = sl_ex_coeff_var_update(sl_block, new_exvar, new_coeff, new_obsVar, new_transVar)
%%Function Description
% This function updates the fields mu, coeff, Sigma, and Omega for an object of
% the class SLBlock (see SLBlock.m for details)

sl_block.exvar = new_exvar;
sl_block.coeff = new_coeff;
sl_block.obsVar = new_obsVar;
sl_block.transVar = new_transVar;

sl_block.H_coeff = sl_make_Hmat(sl_block.T, sl_block.coeff.trans_lag, sl_block.coeff.trans_lead, 1);
sl_block.G_coeff = sl_make_Hmat(sl_block.T, -1*sl_block.coeff.obs_lag, -1*sl_block.coeff.obs_lead, sl_block.coeff.zero);

sl_block.X_tilde = sl_make_Xtilde(sl_block.T, sl_block.coeff.trans_lag, sl_block.coeff.trans_lead);
sl_block.X_tilde_obs = sl_make_Xtilde(sl_block.T, sl_block.coeff.obs_lag, sl_block.coeff.obs_lead);
sl_block.alpha = sl_block.X_tilde * sl_block.ic;
sl_block.xi = sl_block.X_tilde_obs * sl_block.ic_obs;

sl_block.X = sl_block.H_coeff \ sl_block.X_tilde;

end