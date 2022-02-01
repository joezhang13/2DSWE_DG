%Calculate the residual and compute the time derivative -M^(-1)*R
function [F, Force] = DuDt(U, iM, resdata)

%Constants and the input data
g = 9.8;
rho = 1000;
node = resdata.node;
E2N = resdata.E2N;
I2E = resdata.I2E;
B2E = resdata.B2E;
In = resdata.In;
Bn = resdata.Bn;
nBuilding = resdata.nBuilding;
ne = resdata.ne;
np = resdata.np;
wq1d = resdata.wq1d;
wq2d = resdata.wq2d;
phi = resdata.phi;
dphidxi = resdata.dphidxi;
phif = resdata.phif;
ni = size(I2E, 1);
nb = size(B2E, 1);
nq = length(wq1d);
R = zeros(ne*np, 3);
Force = zeros(1, 2, nBuilding);

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
    [flux, ~] = RoeFlux(uqL', flip(uqR, 1)', repmat(normal', 1, nq));   %clockwise for the R element
    Fhat = flux';
    ip = E2N(elemL, :); ip(faceL) = [];
    p1 = node(ip(1), :); p2 = node(ip(2), :);
    dl = norm(p1 - p2);                        %edge length
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
    hq = uq(:, 1);
    Fb = zeros(nq, 3);                         %wall boundary flux
    Fb(:, 2) = 0.5*g*hq.^2*normal(1);
    Fb(:, 3) = 0.5*g*hq.^2*normal(2);
    ip = E2N(elem, :); ip(face) = [];
    p1 = node(ip(1), :); p2 = node(ip(2), :);
    dl = norm(p1 - p2);                        %edge length
    R(indice, :) = R(indice, :) + phif(:, :, face)'*diag(wq1d)*Fb*dl;
    %Calculate the force on each building
    bgroup = B2E(i, 3);
    if bgroup > 1
        iB = bgroup - 1;         %index of the building(based on the order in the .gri file)
        FiB = wq1d*Fb(:, 2:3)*rho*dl;
        Force(1, :, iB) = Force(1, :, iB) + FiB;
    end
end

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