function X_phi = ar_make_Xphi(s, ic)
%%Function Description
% Creates a sparse T x p \tilde{X} matrix to sample initial values of a
% state variable block.
%
% Inputs
%        s:   T x 1 vector of state variable
%        ic:  p x1 vector of initial conditions
%
% Output:
%        X_tilde: required T x p matrix

p = max(size(ic));

X_phi = lagmatrix(s, 1:p);

for i = 1:p
   X_phi(1:i, i) = flip(ic(1:i));
end

end