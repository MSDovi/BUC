function X_tilde = df_make_Xtilde(T, Phi)
%%Function Description
% Creates a sparse Tq x pq \tilde{X} matrix to sample initial values of a
% state variable block.
%
% Inputs
%        T:   number of observations
%        Phi: a q x q x p array of coefficients 
%             (see SBlock for documentation)
% Output:
%        X_tilde: required sparse Tq x pq matrix

[junk, q, p] = size(Phi);

X_tilde = zeros(T*q, p*q);
Phi = reshape(permute(Phi, [2 1 3]), size(Phi, 2), [])';

for i = 1:p
    X_tilde(1 : (q*(p-i+1)), (q*(i-1) + 1) : q*i) = Phi( (q*(i-1) + 1) : (p*q), : );
end

X_tilde = sparse(X_tilde);

end