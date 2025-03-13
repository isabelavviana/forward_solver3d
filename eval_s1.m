function [u] = eval_s1(S, zk, Q, sig, tol)

    srcuse = [];
    srcuse.sources = S.r(:,:);
    srcuse.charges = (sig(:).*S.wts(:)).'/4/pi;
    
    pg = 2;
    U = hfmm3d(tol, zk, srcuse, pg);
 
    u = U.pot; 
    u = u(:);
    % this calcualtes S and fix the quadrature.
    

    ixyzs = S.ixyzs(:);
    npatches = S.npatches;

    istarts = Q.row_ptr(1:end-1);
    iends = Q.row_ptr(2:end)-1;
    isrcinds = cell(npatches,1);
    for i=1:npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
    end

    for i=1:S.npts
        iinds = horzcat(isrcinds{Q.col_ind(istarts(i):iends(i))});
        srcinfo = [];
        srcinfo.r = S.r(:,iinds);
        srcinfo.n = S.n(:,iinds);
        targinfo = [];
        targinfo.r = S.r(:,i);
        targinfo.n = S.n(:,i);
        submat = helm3d.kern(zk, srcinfo, targinfo, 's');
        submat(isnan(submat)) = 0;
        submat(isinf(submat)) = 0;
        u(i) = u(i) - submat*(srcuse.charges(iinds).')*4*pi;
    end

    u = u + Q.spmat*sig;

end