function [u_s, partialu] = fwd_solver(tol, zk, d, sensors, mesh)


d = d/norm(d);  % inc direction (unit norm vec)
eta = 0.5;

S = mesh;

[uinc,graduinc] = helm3d.planewave(zk,d,S.r); 

Q = helm3d.neumann.get_quadrature_correction(S, tol, zk, 0.0, S, 'rpcomb-bc');
Q_s = helm3d.dirichlet.get_quadrature_correction(S, tol, zk, [1,0]);

Q2 = Q;
Q2.wnear = Q.wnear(1,:).';
Q2.spmat = conv_rsc_to_spmat(S, Q2.row_ptr, Q2.col_ind, Q2.wnear);

Q2_s = Q_s;
Q2_s.wnear = Q_s.wnear;
Q2_s.spmat = conv_rsc_to_spmat(S, Q2_s.row_ptr, Q2_s.col_ind, Q2_s.wnear);

b = (graduinc(1,:).*S.n(1,:) + graduinc(2,:).*S.n(2,:) + graduinc(3,:).*S.n(3,:)) - (1i*eta*(uinc.'));
b = b.';

A = @(sig)eval_greeneq(S,zk,Q2_s,Q2,tol,eta,sig);
partialu = gmres(A,b);

targ = [10 0 0];
pg = 1;
pgt = 1;
srcuse = [];
srcuse.sources = S.r(:,:);
srcuse.charges = (partialu(:).*S.wts(:)).'/4/pi;
u_s = hfmm3d(tol, zk, srcuse, pg, targ.', pgt);

end

%to test:
%[u_s, partialu] = fwd_solver(1e-5, 1.1, [1,4,-1], [[10,0,0].',[20,0,0].'], geometries.sphere(1,4))