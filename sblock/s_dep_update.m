function s_block = s_dep_update(s_block, new_dep)
%%Function Description
%
% This function updates y_t for 'objects' of the SBlock 'class' (see
% SVBlock.m for a description of this 'class')
% Inputs:
%   s_block: SBlock 'object', see SBlock.m for a description.
%   new_dep:  T x 1 vector of new dependent vairable; replaces the old
%             dependent vairable in s_block
% Output: 
%   s_block: SBlock 'object' with updated dependent variable

%% Update dependent variable

if size(new_dep) == size(s_block.dep)
    s_block.dep = new_dep;
else

    error('SBlock: Dimensions of old and updated variable do not match.')

end