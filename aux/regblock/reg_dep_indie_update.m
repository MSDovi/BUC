function reg_block = reg_dep_indie_update(reg_block, new_dep, new_indie)
%%Function Description
%
% This function updates y_t and ic for 'objects' of the RegBlock 'class' 
% (see RegBlock.m for a description of this 'class')
% Inputs:
%   reg_block: RegBlock 'object', see RegBlock.m for a description.
%   new_dep:   T x 1 vector of new dependent variable; replaces the old
%              dependent vairable in reg_block
%   new_indie: T x k matrix of new independent variable replaces the old
%              independent vairable in reg_block
% Output: 
%   reg_block: ARBlock 'object' with updated dependent variable

%% Update dependent and independent variable
reg_block.dep = new_dep;
reg_block.indie = new_indie; 
    
end