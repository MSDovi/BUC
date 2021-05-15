function ist = is_stat(phis)
%%Function Description
%
% This function checks the stationarity of AR process with coefficients phis
% 
% Input: 
%   phis: p x 1 vector of autoregressive coefficients
% Output: 
%   ist: Boolean whether AR is stationary

p = max(size(phis));
F = [eye(p-1) zeros(p-1, 1)];
F = [phis'; F];
ist = all(abs(eig(F)) < 1);

end