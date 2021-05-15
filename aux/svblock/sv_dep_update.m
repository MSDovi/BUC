function sv_block = sv_dep_update(sv_block, new_dep)
%%Function Description
%
% This function updates y_t for 'objects' of the SVBlock 'class' (see
% SVBlock.m for a description of this 'class')
% Inputs:
%   sv_block: SVBlock 'object', see SVBlock.m for a description.
%   new_dep:  T x 1 vector of new dependent vairable; replaces the old
%             dependent vairable in sv_block
% Output: 
%   sv_block: SVBlock 'object' with updated dependent variable

%% Update dependent variable

if size(new_dep) == size(sv_block.dep)
    sv_block.dep = new_dep;
else

    error('SVBlock: Dimensions of old and updated variable do not match.')

end