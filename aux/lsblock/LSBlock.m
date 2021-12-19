function ls_block = LSBlock(dep_init, opt)

[T jk] = size(dep_init);
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
HH = H'*H;
HH = H'*H;
H2 = H*H;
H2H2 = H2'*H2;
Pgam = [HH sparse(T,T); sparse(T,T) H2H2];

iVbeta = sparse(zeros(4, 4));
iVbeta(1, 1) = opt.pVic_trend;
iVbeta(2, 2) = opt.pVic_drift;
iVbeta(3, 3) = opt.pVar_sig_trend;
iVbeta(4, 4) = opt.pVar_sig_drift;

var_init = 1/(opt.S * (opt.nu - 1));

ls_block.nu_sig0 = opt.nu;
ls_block.S = opt.S;
ls_block.tau0 = opt.pic_trend;
ls_block.alp0 = opt.pic_drift;
ls_block.sig_trend = opt.pmean_sig_trend;
ls_block.sig_drift = opt.pmean_sig_drift;
ls_block.Pgam = Pgam;
ls_block.H = H;
ls_block.dep = dep_init;
ls_block.iVbeta = iVbeta;
ls_block.beta0 = [opt.pic_trend; opt.pic_drift; opt.pmean_sig_trend; opt.pmean_sig_drift];
ls_block.sig_obs = var_init;
ls_block.T = T;




end