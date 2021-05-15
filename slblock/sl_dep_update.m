function sl_block = sl_dep_update(sl_block, new_dep)
%%Function Description
%
% This function updates y_t for 'objects' of the SLBlock 'class' (see
% SLBlock.m for a description of this 'class')
% Inputs:
%   sl_block: SLBlock 'object', see SLBlock.m for a description.
%   new_dep:  T x 1 vector of new dependent vairable; replaces the old
%             dependent vairable in s_block
% Output: 
%   sl_block: SLBlock 'object' with updated dependent variable

%% Update dependent variable

if size(new_dep) == size(sl_block.dep)
    sl_block.dep = new_dep;
else

    error('SLBlock: Dimensions of old and updated variable do not match.')

end