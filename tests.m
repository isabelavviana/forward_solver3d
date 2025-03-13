clear all
clc
% code with all the forward tests look at this one here. More up-to-date.

r = 1;
na = 4;

zk = 1.1;

tol = 1e-7;
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

eta = 0.5;

%%
%%%%%%%%%%%%%% Setting integral operators %%%%%%%%%%%%%% 
% A = I/2+Sp-ietaS
% A1 = S
% A2 = I/2+Sp
A = @(sig)eval_greeneq(S,zk,Q2_s,Q2,tol,eta,sig);
A1 = @(sig)eval_s1(S, zk, Q2_s, sig, tol);
A2 = @(sig)eval_sprime1(S, zk, Q2, sig, tol);
A21 = @(sig)(eval_sprime1(S, zk, Q2, sig, tol)-sig); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%% Testing operator S %%%%%%%%%%%%%%
fprintf('Testing S with inside source\n')
src = [0 0 0];
targ = [10 0 0];
% exact field is exp(ik|targ-src|)/(4pi|targ-src|)
ts = norm(targ-src);
uex = exp(1i*zk*ts)/(4*pi*ts);

tbd = sqrt(S.r(1,:).^2+S.r(2,:).^2+S.r(3,:).^2);
ubd = exp(1i*zk*tbd)./(4*pi*tbd);
ubd = ubd.';

densu = gmres(A1,ubd,20,1e-6,100);
pg = 1;
pgt = 1;
srcuse = [];
srcuse.sources = S.r(:,:);
srcuse.charges = (densu(:).*S.wts(:)).'/4/pi;
valtarg = hfmm3d(tol, zk, srcuse, pg, targ.', pgt);

fprintf('Error of operator S=%d\n',norm(valtarg.pottarg-uex)/norm(uex))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%% Test for incoming wave - checking simmetry %%%%%%%%%%%%%%
fprintf('Test for incoming wave - checking simmetry and if the results are the same for all the operators\n')
dir = [-1,4,1]; dir = dir/norm(dir);  % inc direction (unit norm vec)
[uinc,graduinc] = helm3d.planewave(zk,dir,S.r); 
b = (graduinc(1,:).*S.n(1,:) + graduinc(2,:).*S.n(2,:) + graduinc(3,:).*S.n(3,:)) - (1i*eta*(uinc.'));
b = b.';
b1 = uinc;
b2 = (graduinc(1,:).*S.n(1,:) + graduinc(2,:).*S.n(2,:) + graduinc(3,:).*S.n(3,:)).';


partialu = gmres(A,b,20,1e-6,100);
partialu1 = gmres(A1,b1,20,1e-6,100);
partialu2 = gmres(A2,b2,20,1e-6,100);
fprintf('Error S=%d\n',norm(partialu1- partialu)./norm(partialu))
fprintf('Error Sp=%d\n',norm(partialu2- partialu)./norm(partialu))

%checking with Manas solver
rhs = -uinc;
pt.r = [10;0;0];
new_dens = helm3d.solver(S, 'dir', rhs, 1e-06, zk);  
new_upt = helm3d.eval(S, 'dir', new_dens, pt, 1e-06, zk);       % eval u scatt @ pt
srcuse.charges = (partialu1(:).*S.wts(:)).'/4/pi;
U_targ1_new = hfmm3d(tol, zk, srcuse, pg, targ.', pgt);
fprintf('Error manas=%d\n',norm(new_upt - (-U_targ1_new.pottarg))./norm(new_upt))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%% Testing operator S' %%%%%%%%%%%%%%
fprintf('Testing Sprime with inside source\n')
src = [0 0 0];
targ = [10 0 0];
% exact field is exp(ik|targ-src|)/(4pi|targ-src|)
ts = norm(targ-src);
uex = exp(1i*zk*ts)/(4*pi*ts);

tbd = sqrt(S.r(1,:).^2+S.r(2,:).^2+S.r(3,:).^2);
ubd = exp(1i*zk*tbd)./(4*pi*tbd);

ubd = ubd.';

% questions:
%  -I/2 + S'(partial u) = partial uinc
% to find densu we need to do gmres (A2, partial uinc + I/2) instead of
% (A2,ubd)
% So ubd should be partial uinc and then we add sigma/2?

% partial uinc is 

densu2 = gmres(A21,ubd,20,1e-6,100); 
pg = 1;
pgt = 1;
srcuse = [];
srcuse.sources = S.r(:,:);
srcuse.charges = (densu2(:).*S.wts(:)).'/4/pi;
valtarg2 = hfmm3d(tol, zk, srcuse, pg, targ.', pgt);
%

fprintf('Error of operator Sprime=%d\n',norm(valtarg.pottarg-uex)/norm(uex))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% evaluate field at target points using green's eq 
% Su=ui, derivative of this equation and combination

pg = 1;
pgt = 1;
targ = [10 10; 10 -10;0 0];
srcuse = [];
srcuse.sources = S.r(:,:);
srcuse.charges = (partialu(:).*S.wts(:)).'/4/pi;
U_targ = hfmm3d(tol, zk, srcuse, pg, targ, pgt);

srcuse.charges = (partialu1(:).*S.wts(:)).'/4/pi;
U_targ1 = hfmm3d(tol, zk, srcuse, pg, targ, pgt);


srcuse.charges = (partialu2(:).*S.wts(:)).'/4/pi;
U_targ2 = hfmm3d(tol, zk, srcuse, pg, targ, pgt);


% u_s = eval_s_target(S, zk, Q2, partialu, U_targ, srcuse, targ);


