function ls_block = ls_dep_update(ls_block, new_dep)
%%Function Description
%
% This function updates y_t and ic for 'objects' of the LSBlock 'class' 
% (see LSBlock.m for a description of this 'class')
% Inputs:
%   ls_block:  LSBlock 'object', see LSBlock.m for a description.
%   new_dep:   T x 1 vector of new dependent variable; replaces the old
%              dependent vairable in ls_block
% Output: 
%   ls_block:  LSBlock 'object' with updated dependent variable

%% Update dependent 
ls_block.dep = new_dep;
    
end