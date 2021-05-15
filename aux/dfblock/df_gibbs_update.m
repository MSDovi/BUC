function df_block = df_gibbs_update(df_block)
%%Function Description
%
% This function samples f_t in the model described in DFBlock.m
% 
% Input: 
%   df_block: DFBlock 'object', see DFBlock.m for a description.
% Output: 
%   df_block: DFblock 'object' with updated f, f0, alpha

%% Setup
y = df_block.dep;
T = df_block.T;
iSigma = inv(df_block.obsVar);
iOmega = inv(df_block.transVar);
H = df_block.H_coeff;
c0 = H\(df_block.alpha + df_block.exvar);
X = df_block.X;
iV_0 = inv(df_block.opt.icV);
ic_0 = df_block.opt.ic0;
p = df_block.p;
Lambda = df_block.loadings;
Lambda_s = df_block.loadings_s;
q = df_block.q;
n = df_block.n;
VLambda = df_block.opt.Vl;
MLambda = df_block.opt.ml;
nu = df_block.opt.nu;
S = df_block.opt.S;

%% sample f
HiOgH = H'*iOmega*H;
Kf = HiOgH + Lambda'*iSigma*Lambda;

f_hat = Kf\(HiOgH*c0 + Lambda'*iSigma*y);
f = f_hat + chol(Kf,'lower')'\randn(T*q,1);
F = reshape(f, q, T)';
%% sample f0
f_tilde = f - H\df_block.exvar;
Kf0 = X'*HiOgH*X + iV_0;
fhat = Kf0\(X'*HiOgH*f_tilde + iV_0*ic_0);
f0 = fhat + chol(Kf0,'lower')'\randn(p*q,1);

%% sample Lambda -- equation by equation
count_Lam = 0;
for isave = 2:n 
    slct = isave:n:((T-1)*n+isave);
    Sigma_j = df_block.obsVar(slct, slct);
    y_j = y(slct);
    if isave <= q
        q_i = isave-1;
        sp = count_Lam+1:count_Lam+q_i;
        Xf = F(:,1:q_i); 
        K_Lami = inv(VLambda(sp, sp)) ...
            + Xf'*(Sigma_j\Xf);
        Lami_hat = K_Lami\(inv(VLambda(sp, sp))*MLambda(sp) + ...
                    Xf'*(Sigma_j\(y_j-F(:,isave)))); 
    else
        q_i = q;
        sp = count_Lam+1:count_Lam+q_i;
        Xf = F;
        K_Lami = inv(VLambda(sp, sp)) ...
            + Xf'*(Sigma_j\Xf);
        Lami_hat = K_Lami\(inv(VLambda(sp, sp))*MLambda(sp) + ... 
                            Xf'*(Sigma_j\y_j));            
    end
    Lambda_s(isave,1:q_i) = Lami_hat + chol(K_Lami,'lower')'\randn(q_i,1);
    count_Lam = count_Lam + q_i;
end    

%% Save some preliminary data
df_block.f = f;
df_block.F = F;
df_block.ic = f0;
df_block.F0 = reshape(df_block.ic, df_block.p, df_block.q);
df_block.alpha = df_block.X_tilde * df_block.ic;
df_block.res = df_block.H_coeff*df_block.f - df_block.alpha - df_block.exvar;
df_block.Res = reshape(df_block.res, df_block.q, df_block.T)';
df_block.Res_obs = reshape(y - Lambda*f, df_block.n, df_block.T)';
df_block.loadings = kron(speye(df_block.T), Lambda_s);
df_block.loadings_s = Lambda_s;

%% Sample homoscedastic variance for observation equation
Sig_y = zeros(n, 1);
for i = 1:n
    newS_y = S +  sum(df_block.Res_obs(:, i).^2)'/2;
    Sig_y(i) = 1./gamrnd(nu + T/2, 1./newS_y);    
end

%% Save homoscedastic variance

df_block.obsVar = kron(speye(df_block.T), sparse(diag(Sig_y)));

end