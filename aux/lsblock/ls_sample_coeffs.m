function [beta beta_hat Kbeta] = sample_beta(y,X,sigy2,beta0,iVbeta)
Kbeta = iVbeta + X'*X/sigy2 ; 
beta_hat = Kbeta\(iVbeta*beta0 + X'*y/sigy2);
beta = beta_hat + chol(Kbeta,'lower')'\randn(4,1);
end
