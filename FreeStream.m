%Free stream and free-stream preservation test
function [normR, dt] = FreeStream(mesh, CFL0, Nstep, p)

%Mesh information
node = mesh.Node;
E2N = mesh.Elem;
ne = mesh.nElem;
[I2E, B2E, In, Bn, ~] = processMesh(mesh);

%Free-stream initial condition and time step
np = (p + 1)*(p + 2)/2;
hf = 1.25; huf = 4.5; hvf = -6.7;
% hf = 1.2; huf = 0.1; hvf = 0.2;
uf = [hf, huf, hvf];
U = repmat(uf, ne*np, 1);
Dt = zeros(ne, 1);      %time step for each cell
Rot = [0, 1; -1, 0];      %rotation matrix
s = zeros(3, 1);
for i = 1 : ne
    pp = node(E2N(i, :), :);
    l1 = pp(3, :) - pp(2, :); n1 = Rot*l1'; n1 = n1/norm(n1);
    l2 = pp(1, :) - pp(3, :); n2 = Rot*l2'; n2 = n2/norm(n2);
    l3 = pp(2, :) - pp(1, :); n3 = Rot*l3'; n3 = n3/norm(n3);
    P = norm(l1) + norm(l2) + norm(l3);
    A = 0.5*abs(l1(1)*l2(2) - l2(1)*l1(2));
    d = 2*A/P;
    [~, s(1)] = RoeFlux(uf', uf', n1);
    [~, s(2)] = RoeFlux(uf', uf', n2);
    [~, s(3)] = RoeFlux(uf', uf', n3);
    smax = max(s);
    Dt(i) = CFL0/(p + 1)*d/smax;
end
dt = min(Dt);

%Pre-computed quantities
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
resdata.ne = ne;
resdata.np = np;
resdata.wq1d = wq1d;
resdata.wq2d = wq2d;
resdata.phi = phi;
resdata.dphidxi = dphidxi;
resdata.phif = phif;
resdata.uf = uf;

%Loop over time steps
normR = zeros(Nstep, 1);
for n = 1 : Nstep
    %Update the states
    [F0, R] = DuDtFS(U, iM, resdata);
    [F1, ~] = DuDtFS(U + 0.5*dt*F0, iM, resdata);
    [F2, ~] = DuDtFS(U + 0.5*dt*F1, iM, resdata);
    [F3, ~] = DuDtFS(U + dt*F2, iM, resdata);
    U = U + dt/6*(F0 + 2*F1 + 2*F2 + F3);
    normR(n) = R;
    disp(['Step ', num2str(n), '. L1 norm of the residual: ', num2str(normR(n))]);
end

end

%Calculate the residual (and L1 norm) and compute the time derivative -M^(-1)*R
function [F, normR] = DuDtFS(U, iM, resdata)

%Constants and the input data
g = 9.8;
node = resdata.node;
E2N = resdata.E2N;
I2E = resdata.I2E;
B2E = resdata.B2E;
In = resdata.In;
Bn = resdata.Bn;
ne = resdata.ne;
np = resdata.np;
wq1d = resdata.wq1d;
wq2d = resdata.wq2d;
phi = resdata.phi;
dphidxi = resdata.dphidxi;
phif = resdata.phif;
uf = resdata.uf;
ni = size(I2E, 1);
nb = size(B2E, 1);
nq = length(wq1d);
R = zeros(ne*np, 3);

%Loop over element interiors
for i = 1 : ne
    pp = node(E2N(i, :), :);
    x1 = pp(1, :)'; x2 = pp(2, :)'; x3 = pp(3, :)';
    J = [x2(1)-x1(1), x3(1)-x1(1); x2(2)-x1(2), x3(2)-x1(2)];          %mapping Jacobian
    iJ = [x3(2)-x1(2), x1(1)-x3(1); x1(2)-x2(2), x2(1)-x1(1)]/det(J);  %inverse of mapping Jacobian
    indice = (i-1)*np+1 : i*np;
    u = U(indice, :);
    uq = phi*u;             %states at the quadrature points
    uq1 = uq(:, 1); uq2 = uq(:, 2); uq3 = uq(:, 3);
    Fx = [uq2, uq2.^2./uq1+0.5*g*uq1.^2, uq2.*uq3./uq1];            %fluxes at the quadrature points
    Fy = [uq3, uq2.*uq3./uq1, uq3.^2./uq1+0.5*g*uq1.^2];            %fluxes at the quadrature points
    gradx = dphidxi(:, :, 1)*iJ(1, 1) + dphidxi(:, :, 2)*iJ(2, 1);  %gradient matrix in global space
    grady = dphidxi(:, :, 1)*iJ(1, 2) + dphidxi(:, :, 2)*iJ(2, 2);  %gradient matrix in global space
%     temp = zeros(ne*np, 3);
%     temp(indice, :) = -(gradx'*diag(wq2d)*Fx + grady'*diag(wq2d)*Fy)*det(J);
%     R = R + temp;
    R(indice, :) = R(indice, :) - (gradx'*diag(wq2d)*Fx + grady'*diag(wq2d)*Fy)*det(J);
end

%Loop over interior faces
for i = 1 : ni
    elemL = I2E(i, 1); faceL = I2E(i, 2); 
    elemR = I2E(i, 3); faceR = I2E(i, 4);
    normal = In(i, :);
    indiceL = (elemL-1)*np+1 : elemL*np;
    indiceR = (elemR-1)*np+1 : elemR*np;
    uL = U(indiceL, :); uR = U(indiceR, :);
    uqL = phif(:, :, faceL)*uL;
    uqR = phif(:, :, faceR)*uR;
    
    [flux, ~] = RoeFlux(uqL', flip(uqR, 1)', repmat(normal', 1, nq));  %clockwise for the R element
    Fhat = flux';
%     Fhat = zeros(nq, 3);
%     for j = 1 : nq
%         [flux, ~] = RoeFlux(uqL(j, :)', uqR(nq+1-j, :)', normal');  %clockwise for the R element
%         Fhat(j, :) = flux';
%     end

    ip = E2N(elemL, :); ip(faceL) = [];
    p1 = node(ip(1), :); p2 = node(ip(2), :);
    dl = norm(p1 - p2);                        %edge length
    
%     temp = zeros(ne*np, 3);
%     temp(indiceL, :) = temp(indiceL, :) + phif(:, :, faceL)'*diag(wq1d)*Fhat*dl;
%     temp(indiceR, :) = temp(indiceR, :) - phif(:, :, faceR)'*diag(wq1d)*flip(Fhat, 1)*dl;
%     R = R + temp;
    R(indiceL, :) = R(indiceL, :) + phif(:, :, faceL)'*diag(wq1d)*Fhat*dl;
    R(indiceR, :) = R(indiceR, :) - phif(:, :, faceR)'*diag(wq1d)*flip(Fhat, 1)*dl;
end

%Loop over boundary faces
for i = 1 : nb
    elem = B2E(i, 1); face = B2E(i, 2);
    normal = Bn(i, :);
    indice = (elem-1)*np+1 : elem*np;
    u = U(indice, :);
    uq = phif(:, :, face)*u;
    
    [flux, ~] = RoeFlux(uq', repmat(uf', 1, nq), repmat(normal', 1, nq)); %free-stream boundary condition
    Fb = flux';
%     Fb = zeros(nq, 3);                         %free-stream boundary condition
%     for j = 1 : nq
%         [flux, ~] = RoeFlux(uq(j, :)', uf', normal');
%         Fb(j, :) = flux';
%     end
    
    ip = E2N(elem, :); ip(face) = [];
    p1 = node(ip(1), :); p2 = node(ip(2), :);
    dl = norm(p1 - p2);                        %edge length
    
%     temp = zeros(ne*np, 3);
%     temp(indice, :) = phif(:, :, face)'*diag(wq1d)*Fb*dl;
%     R = R + temp;
    R(indice, :) = R(indice, :) + phif(:, :, face)'*diag(wq1d)*Fb*dl;
end

%Calculate the residual norm
Res = reshape(R', 3*ne*np, 1);
normR = norm(Res, 1);

%Evaluate F = -M^(-1)*R
F = zeros(ne*np, 3);
for i = 1 : ne
    pp = node(E2N(i, :), :);
    x1 = pp(1, :)'; x2 = pp(2, :)'; x3 = pp(3, :)';
    J = [x2(1)-x1(1), x3(1)-x1(1); x2(2)-x1(2), x3(2)-x1(2)];
    indice = (i-1)*np+1 : i*np;
    F(indice, :) = -iM/det(J)*R(indice, :);
end

end