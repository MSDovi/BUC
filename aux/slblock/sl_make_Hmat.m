function H = sl_make_Hmat(T, a, b, c)
%%Function Description
% Creates a sparse T x T matrix with c_i on the diagonal, and -a_i on the ith
% lower off-diagonal, and -b_i on the ith upper off-diagonal
%
% Inputs
%        T: size of matrix
%        c: T x 1 vector or scalar to go on diagonal
%        a: g x 1 vector
%            a(0) contains the value to be put on the first lower-diagonal
%            a(1) contains the value to be put on the second lower-diagonal
%            a(g-1) contains the value to be put on the gth lower-diagonal
%        b: analogous to a
% Output:
%        H: required sparse T x T matrix

g = max(size(a));
h = max(size(b));
if g >= T
    error('Not enough observations to make H matrix.')
end
H = sparse(1:T, 1:T, c);

for i=1:g
    H = H - a(i)*sparse((i+1):T, 1:(T-i), ones(1, T-i), T, T);
end

for i = 1:h
    H = H - (b(i)*sparse((i+1):T, 1:(T-i), ones(1, T-i), T, T))';
end
end