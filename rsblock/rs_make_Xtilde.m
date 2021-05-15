function X_tilde = rs_make_Xtilde(T, phi)
%%Function Description
% Creates a sparse T x p \tilde{X} matrix to sample initial values of a
% state variable block.
%
% Inputs:
%        T:   number of observations
%        phi: a p x 1 vector of coefficients (see SBlock for documentation)
% Output:
%        X_tilde: required sparse T x p matrix

p = max(size(phi));

X_tilde = sparse(zeros(T, p));

for i = 1:p
   X_tilde(:, i) = [phi(i:p)' zeros(1, T - p + i - 1)]';
end

X_tilde = sparse(X_tilde);

end