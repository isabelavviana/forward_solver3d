clear all
clc

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
jn = @(n,z) sqrt(pi/2/z)*besselj(n+0.5,z);
hn = @(n,z) sqrt(pi/2/z)*besselh(n+0.5,1,z);

jnp = @(n,z) 0.5*(jn(n-1,z) - (jn(n,z) + z*jn(n+1,z))/z);
hnp = @(n,z) 0.5*(hn(n-1,z) - (hn(n,z) + z*hn(n+1,z))/z);

zfac = 1j*zk*zk*(jn(1,zk)*hnp(1,zk) + jnp(1,zk)*hn(1,zk))/2;
zfac2 = 1j*zk*jn(1,zk)*hn(1,zk);
%%
rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
sig = S.r(3,:).'./rr.';
eta = 0.5;

% u = eval_sprime(S, zk, Q2, sig, tol,U);
% u2 = eval_s(S,zk,Q2_s,sig,tol,U);
[u_tot,u2,u] = eval_greeneq(S,zk,Q2_s,Q2,sig,tol,eta);



%%
uex = sig.*zfac;
err = norm((u-uex).*sqrt(S.wts(:)))/norm(uex.*sqrt(S.wts(:)));
fprintf('error in sprime=%d\n',err);

%%
uex2 = sig.*zfac2;
err2 = norm((u2-uex2).*sqrt(S.wts(:)))/norm(uex2.*sqrt(S.wts(:)));
fprintf('error in s=%d\n',err2);

%%
% u_tot_ex = sig/2 + uex - (sqrt(-1)*eta*uex2);
% err = norm((u_tot-u_tot_ex).*sqrt(S.wts(:)))/norm(u_tot_ex.*sqrt(S.wts(:)));
% fprintf('error in total u=%d\n',err);

function [u] = eval_sprime(S, zk, Q, sig, tol,U,srcuse)
    
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

    u = u + Q.spmat*sig;

end



function [u] = eval_s(S, zk, Q, sig, tol,U,srcuse)
 
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
    u = u + Q.spmat*sig;

end



function [u_tot,u_s,u_sprime] = eval_greeneq(S,zk,Q1,Q2,sig,tol,eta,U)

srcuse = [];
srcuse.sources = S.r(:,:);
srcuse.charges = (sig(:).*S.wts(:)).'/4/pi;
pg = 2;
U = hfmm3d(tol, zk, srcuse, pg);

    
u_s = eval_s(S,zk,Q1,sig,tol,U,srcuse);
u_sprime = eval_sprime(S, zk, Q2, sig, tol,U,srcuse);

u_tot = sig/2 + u_sprime - (sqrt(-1)*eta*u_s);

end