function ar_block = ar_dep_ic_update(ar_block, new_dep, new_ic)
%%Function Description
%
%This function updates y_t and ic for 'objects' of the ARBlock 'class' (see
%ARBlock.m for a description of this 'class')
%
%Inputs:
%   ar_block: ARBlock 'object', see ARBlock.m for a description.
%   new_dep:  T x 1 vector of new dependent vairable; replaces the old
%             dependent vairable in ar_block
%   new_ic:   p x 1 vector of new initial conditions
%
%Output: 
%   s_block:  ARBlock 'object' with updated dependent variable

%% Update dependent variable
if size(new_dep) == size(ar_block.dep)
    ar_block.dep = new_dep;
else

    error('ARBlock: Dimensions of old and updated variable do not match.')
end

if size(new_ic) == size(ar_block.ic)
    ar_block.ic = new_ic;
    ar_block.X_phi = ar_make_Xphi(ar_block.dep, ar_block.ic);

else

    error('ARBlock: Dimensions of old and updated variable do not match.')
    
end