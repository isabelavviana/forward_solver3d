function [u] = eval_s_target(S, zk, Q, dudn, U, srcuse, targ)

    u = U.pottarg;  
    u = u(:);  

    ixyzs = S.ixyzs(:);
    npatches = S.npatches;

    istarts = Q.row_ptr(1:end-1);
    iends = Q.row_ptr(2:end)-1;
    isrcinds = cell(npatches,1);

    for i=1:npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
    end

    
    ntarg = size(targ, 2); 
    for j=1:ntarg
        targinfo = [];
        targinfo.r = targ(:,j);  
        targinfo.n = [0; 0; 0]; 

        for i=1:S.npts
            iinds = horzcat(isrcinds{Q.col_ind(istarts(i):iends(i))});
            srcinfo = [];
            srcinfo.r = S.r(:,iinds);
            srcinfo.n = S.n(:,iinds);

            submat = helm3d.kern(zk, srcinfo, targinfo, 's');
            submat(isnan(submat)) = 0;
            submat(isinf(submat)) = 0;

         
            u(j) = u(j) - submat*(srcuse.charges(iinds).')*4*pi;
        end
    end

    u = u + Q.spmat*dudn;

end