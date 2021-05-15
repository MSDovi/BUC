function H = df_make_Hmat(T, Phi)
%%Function Description
% Creates a sparse T x T matrix with 1 on the diagonal, and -a_i on the ith
% lower off-diagonal.
%
% Inputs
%        T: time dimension
%        a: q x q x p array
%            a(:, :, 1) contains matrix for first lower-diagonal
%            a(:, :, 2) contains matrix for second lower-diagonal
%            a(:, :, g-1) contains matrix for gth lower-diagonal
% Output:
%        H: required sparse Tq x Tq matrix

[junk, q, p] = size(Phi);
if p >= T
    error('Not enough observations to make H matrix.')
end

H = speye(T*q);

for i = 1:p

    hzeros = sparse(q*i, (T-i)*q);
    vzeros = sparse(T * q, q*i);
    aux = kron(speye(T-i),Phi(:, :, i));
    aux = cat(1,hzeros,aux);
    H = H - cat(2, aux, vzeros);

end
    
end