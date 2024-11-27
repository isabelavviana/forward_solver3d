%forward_solver Isabela
r = 1;
na = 4;

zk = 1.1;

tol = 1e-5;
S = geometries.sphere(r, na);

Q = helm3d.neumann.get_quadrature_correction(S, tol, zk, 0.0, S, 'rpcomb-bc');
Q_s = helm3d.dirichlet.get_quadrature_correction(S, tol, zk, [1,0]);
%%
Q2 = Q;
Q2.wnear = Q.wnear(1,:).';
Q2.spmat = conv_rsc_to_spmat(S, Q2.row_ptr, Q2.col_ind, Q2.wnear);

Q2_s = Q_s;
Q2_s.wnear = Q_s.wnear;
Q2_s.spmat = conv_rsc_to_spmat(S, Q2_s.row_ptr, Q2_s.col_ind, Q2_s.wnear);
%%
%%

eta = 0.5;

A = eval_greeneq(S,zk,Q2_s,Q2,tol,eta);

dir = [-1,4,1]; dir = dir/norm(dir);  % inc direction (unit norm vec)

uinc = helm3d.planewave(zk,dir,S.r); %b

%AX = b
%GMRES @sigma eval_greeneq(S,zk,Q1,Q2,sig,tol,eta,U) to solve for x


%I'm having trouble with the dimension of A, u_s and u_sprime are vectors
%of 2880x1 but i need matrices operators.

sigma= gmres(A,uinc);

function [u] = eval_sprime(S, zk, Q, tol,U,srcuse)
    
    u = U.grad(1,:).*S.n(1,:) + U.grad(2,:).*S.n(2,:) + ...
        U.grad(3,:).*S.n(3,:);    
    u = u(:);
    %this calcualtes Sprime and fix the quadrature.
    

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

    u = u + Q.spmat;

end



function [u] = eval_s(S, zk, Q, tol,U,srcuse)
 
    u = U.pot; 
    u = u(:);
    % sum(isinf(u))
    %this calcualtes S and fix the quadrature.
    

    ixyzs = S.ixyzs(:);
    npatches = S.npatches;

    istarts = Q.row_ptr(1:end-1);
    iends = Q.row_ptr(2:end)-1;
    isrcinds = cell(npatches,1);
    for i=1:npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
    end

    % iinf = 0;
    % inan = 0;
    for i=1:S.npts
        iinds = horzcat(isrcinds{Q.col_ind(istarts(i):iends(i))});
        srcinfo = [];
        srcinfo.r = S.r(:,iinds);
        srcinfo.n = S.n(:,iinds);
        targinfo = [];
        targinfo.r = S.r(:,i);
        targinfo.n = S.n(:,i);
        submat = helm3d.kern(zk, srcinfo, targinfo, 's');
        % size(submat)
        % inan = inan  + sum(isnan(submat));
        % iinf = iinf  + sum(isinf(submat));
        submat(isnan(submat)) = 0;
        submat(isinf(submat)) = 0;
        u(i) = u(i) - submat*(srcuse.charges(iinds).')*4*pi;
    end
    % inan
    % iinf

    % sum(isinf(u))
    u = u + Q.spmat;

end



function [A] = eval_greeneq(S,zk,Q1,Q2,tol,eta)

rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
sig = S.r(3,:).'./rr.';
srcuse = [];
srcuse.sources = S.r(:,:);
srcuse.charges = (sig(:).*S.wts(:)).'/4/pi;
pg = 2;
U = hfmm3d(tol, zk, srcuse, pg);

    
S_op = eval_s(S,zk,Q1,tol,U,srcuse);
Sprime_op = eval_sprime(S, zk, Q2, tol,U,srcuse);

A = ((eye(size(Sprime_op))/2 - Sprime_op +sqrt(-1)*eta*S_op));

end