function X_tilde = sl_make_Xtilde(T, phi, theta)
%%Function Description
% Creates a sparse T x p \tilde{X} matrix to sample initial values of a
% state variable block.
%
% Inputs
%        T:     number of observations
%        phi:   a q_1 x 1 vector of coefficients (see SBlock for documentation)
%        theta: a q_2 x 1 vector of coefficients (see SBlock for documentation)
% Output:
%        X_tilde: required sparse T x (q_1 + q_2) matrix

q_1 = max(size(phi));
q_2 = max(size(theta));

X_tilde = sparse(zeros(T, q_1));

for i = 1:q_1
   X_tilde(:, i) = [phi(i:q_1)' sparse(1, T - q_1 + i - 1)]';
end

for i = (q_1 + 1):(q_1 + q_2)
   
   X_tilde(:, i) = flip([theta((i - q_1):q_2)' sparse(1, T - q_2 + i - q_1 - 1)]');
end

X_tilde = sparse(X_tilde);

end