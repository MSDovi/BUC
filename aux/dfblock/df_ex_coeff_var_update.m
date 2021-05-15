function df_block = df_ex_coeff_var_update(df_block, new_exvar, new_coeff, new_obsVar, new_transVar)
%%Function Description
% This function updates the fields mu, phi, Sigma, and Omega for an object
% of the class DFBlock (see DFBlock.m for details)

df_block.coeff = new_coeff;
df_block.H_coeff = df_make_Hmat(df_block.T, df_block.coeff);


if size(new_exvar, 1) == df_block.T
    df_block.exvar = reshape(new_exvar',df_block.T*df_block.q, 1);
else
    df_block.exvar = new_exvar;
end

if max(size(new_obsVar)) == df_block.n
    df_block.obsVar = kron(speye(df_block.T), new_obsVar);
elseif max(size(new_obsVar)) == (df_block.n * df_block.T)
    df_block.obsVar = new_obsVar;
end

if max(size(new_transVar)) == df_block.q
    df_block.transVar = kron(speye(df_block.T), new_transVar);
elseif max(size(new_transVar)) == (df_block.q * df_block.T)
    df_block.transVar = new_transVar;
end

df_block.X_tilde = df_make_Xtilde(df_block.T, df_block.coeff);
df_block.X = df_block.H_coeff \ df_block.X_tilde;
df_block.Exvar = reshape(df_block.exvar, df_block.q, df_block.T)';

end