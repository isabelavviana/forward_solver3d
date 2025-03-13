function [out] = eval_greeneq(S,zk,Q1,Q2,tol,eta,sig)

    srcuse = [];
    srcuse.sources = S.r(:,:);
    srcuse.charges = (sig(:).*S.wts(:)).'/4/pi;
    
    pg = 2;
    U = hfmm3d(tol, zk, srcuse, pg);
    
    s_op = eval_s(S,zk,Q1,sig,U,srcuse);
    s_prime = eval_sprime(S, zk, Q2, sig,U,srcuse);
    
    out = sig/2 + s_prime - (1i*eta*s_op);

end
