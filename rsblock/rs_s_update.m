function rs_block = rs_s_update(rs_block, new_s)
%%Function Description
% This function updates the s field of an s_block for a 'residual' state
% variable. See RSBlock.m for documentation.

rs_block.s = new_s;

end