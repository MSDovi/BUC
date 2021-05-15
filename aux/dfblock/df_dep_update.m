function df_block = df_dep_update(df_block, new_dep)
%%Function Description
%
% This function updates the field of the dependent variable for objects of
% the class DFBlock (see DFBlock.m for details)

%% Update dependent variable

if size(new_dep, 1) == df_block.T
    df_block.dep = reshape(new_dep',df_block.T*df_block.n, 1);
else
    df_block.dep = new_dep;
end

df_block.Dep = reshape(df_block.dep, df_block.n, df_block.T)';

end