function H = make_Hmat(T, a)
%%Function Description
% Creates a sparse T x T matrix with 1 on the diagonal, and -a_i on the ith
% lower off-diagonal.
%
% Inputs
%        T: size of matrix
%        a: g x 1 vector
%            a(0) contains the value to be put on the first lower-diagonal
%            a(1) contains the value to be put on the second lower-diagonal
%            a(g-1) contains the value to be put on the gth lower-diagonal
% Output:
%        H: required sparse T x T matrix

g = max(size(a));
if g >= T
    error('Not enough observations to make H matrix.')
end
H = speye(T);

for i=1:g
    H = H - a(i)*sparse((i+1):T, 1:(T-i), ones(1, T-i), T, T);   
end
    
end