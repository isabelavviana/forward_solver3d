function [u] = eval_sprime1(S, zk, Q, sig, tol)
% calculates S'phi with the U inside                    
    srcuse = [];
    srcuse.sources = S.r(:,:);
    srcuse.charges = (sig(:).*S.wts(:)).'/4/pi;
    
    pg = 2;
    U = hfmm3d(tol, zk, srcuse, pg);
    
    u = U.grad(1,:).*S.n(1,:) + U.grad(2,:).*S.n(2,:) + ...
        U.grad(3,:).*S.n(3,:);    
    u = u(:);
   
    %this calculates Sprime
    

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
        submat = helm3d.kern(zk, srcinfo, targinfo, 'sprime');
        submat(isnan(submat)) = 0;
        u(i) = u(i) - submat*(srcuse.charges(iinds).')*4*pi;
    end

    u = u + Q.spmat*sig;
    u = u + sig/2;
   

end
