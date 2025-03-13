clear
clc


% testing the derivative
% F(dB+h) = F(dB) + J(hn) (dB -> r = 1, dB+h -> r=1 + h)
% Jh = w
% RHS = -partialu(h.n) = w


tol = 1e-9;
%plane wave
zk = 1.1;
d = [1,4,-1];
d = d / norm(d);  

%receptor
sensors = [10;10;0];

%geometry
mesh = geometries.sphere(1,6);
M = width(sensors);
S = mesh;

%calculating F(D)
[u_s_initial, partial_u] = fwd_solver(tol, zk, d, sensors.', mesh);
% partial_u = partial_u.';
u_s_initial.pottarg

for ii = 1 : 4

    % update
    h = 0.1^ii;

    %geometry
    meshnew=[];
    meshnew = geometries.sphere(1 + h,6); %meshnew


    %calculating F(D+h)
    [u_s_exact, ~] = fwd_solver(tol, zk, d, sensors.', meshnew);
    u_s_exact.pottarg

    R_l = u_s_initial.pottarg - u_s_exact.pottarg
    % R_l(1)
   
    b = -h*partial_u;


    w = helm3d.solver(S, 'dir', b, tol, zk);
    targinfo.r = sensors;
    % targinfo.patch_id = -1;

    Jdq = helm3d.eval(S, 'dir', w, targinfo, tol, zk);
    Jdq
    imag(u_s_initial.pottarg + Jdq.'- u_s_exact.pottarg)

end
