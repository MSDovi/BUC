function rs_block = rs_ex_coeff_var_update(rs_block, new_exvar, new_coeff, new_Var)
%%Function Description
% This function updates the fields phi and Omega for an object of
% the class RSBlock (see RSBlock.m for details)

rs_block.exvar = new_exvar;
rs_block.coeff = new_coeff;
rs_block.Var = new_Var;

rs_block.X_tilde = rs_make_Xtilde(rs_block.T, rs_block.coeff);
rs_block.X = rs_block.H_coeff \ rs_block.X_tilde;


end