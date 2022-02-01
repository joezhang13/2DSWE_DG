%DG solver for the 2D shallow water equations
function [Uo, Force] = DGSolver(mesh, U0, dt, T, p)

%Mesh information
node = mesh.Node;
E2N = mesh.Elem;
ne = mesh.nElem;
[I2E, B2E, In, Bn, ~] = processMesh(mesh);
nBuilding = mesh.B.nbfgrp - 1;

%Time stepping
Tmax = T(end);
Nt = ceil(Tmax/dt);
dt = Tmax/Nt;

%Pre-computed quantities
np = (p + 1)*(p + 2)/2;
[coeff, ~, ~] = TriLagrange2D(p);
[xq1d, wq1d] = quad1d(2*p + 1);
[xq2d, wq2d] = quad2d(2*p + 1);
phi = basis(xq2d, p, coeff);
dphidxi = gbasis(xq2d, p, coeff);
phif(:, :, 1) = basis([1-xq1d; xq1d], p, coeff);        %basis function matrix for cell faces
phif(:, :, 2) = basis([xq1d*0; flip(xq1d)], p, coeff);
phif(:, :, 3) = basis([xq1d; xq1d*0], p, coeff);
M = phi'*diag(wq2d)*phi;
iM = inv(M);

%Data for residual calculation
resdata.node = node;
resdata.E2N = E2N;
resdata.I2E = I2E;
resdata.B2E = B2E;
resdata.In = In;
resdata.Bn = Bn;
resdata.nBuilding = nBuilding;
resdata.ne = ne;
resdata.np = np;
resdata.wq1d = wq1d;
resdata.wq2d = wq2d;
resdata.phi = phi;
resdata.dphidxi = dphidxi;
resdata.phif = phif;

%Loop over time steps
No = length(T);
Uo = zeros(ne*np, 3, No);
io = 1;
U = U0;
t = 0;
Force = zeros(Nt + 1, 2, nBuilding);          %Forces on buildings(the 3rd dimension is #buildings)
for n = 1 : Nt
    disp(['t = ', num2str(t)]);
    %Output the solution
    if abs(t - T(io)) < dt/2
        Uo(:, :, io) = U;
        if io < No
            io = io + 1;
        end
    end 
    %Update the states
    [F0, Force(n, :, :)] = DuDt(U, iM, resdata);
    [F1, ~] = DuDt(U + 0.5*dt*F0, iM, resdata);
    [F2, ~] = DuDt(U + 0.5*dt*F1, iM, resdata);
    [F3, ~] = DuDt(U + dt*F2, iM, resdata);
    U = U + dt/6*(F0 + 2*F1 + 2*F2 + F3);
    t = t + dt;
end
%Output the final solution
disp(['t = ', num2str(t)]);
Uo(:, :, io) = U;
[~, Force(end, :, :)] = DuDt(U, iM, resdata);

end