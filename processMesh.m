%Generate the matrices of the mesh from a .gri file
function [I2E, B2E, In, Bn, Area] = processMesh(mesh)

%read mesh

%mesh = readgri(fname);
E2N = mesh.Elem;
NB = mesh.B.nodes;
p = mesh.Node';
nBGroup = mesh.B.nbfgrp;
nElemTot = mesh.nElem;
[IE,BE] = edgehash(E2N);

%I2E matrix and In matrix
nIE = size(IE, 1);
I2E = zeros(nIE, 4);
In = zeros(nIE, 2);
for i = 1 : nIE
    if IE(i, 3) < IE(i, 4)
        elemL = IE(i, 3);
        elemR = IE(i, 4);
    else
        elemL = IE(i, 4);
        elemR = IE(i, 3);
    end
    %the local face number equals the local number of the node which is not
    %on the face.
    n1 = IE(i, 1); n2 = IE(i, 2);
    iL = (E2N(elemL, :) ~= n1) & (E2N(elemL, :) ~= n2);
    faceL = find(iL);
    iR = (E2N(elemR, :) ~= n1) & (E2N(elemR, :) ~= n2);
    faceR = find(iR);
    I2E(i, :) = [elemL, faceL, elemR, faceR];
    %find the vector parallel to the cell face and counterclockwise to 
    %elemL, then rotate it by 90 degrees clockwisely
    idx1 = mod(faceL, 3) + 1; idx2 = mod(faceL + 1, 3) + 1;
    faceVec = p(:, E2N(elemL, idx2)) - p(:, E2N(elemL, idx1));
    R = [0, 1; -1, 0];        %rotation matrix
    nVec = R*faceVec;
    nVec = nVec/norm(nVec);   %normal vector pointing from L to R
    In(i, :) = nVec';
end

%B2E matrix and Bn matrix
nBTot = size(BE, 1);
B2E = zeros(nBTot, 3);
Bn = zeros(nBTot, 2);
for i = 1 : nBTot
    elem = BE(i, 3);
    n1 = BE(i, 1); n2 = BE(i, 2);
    iB = (E2N(elem, :) ~= n1) & (E2N(elem, :) ~= n2);
    face = find(iB);
    %determine the boundary group index
    for j = 1 : nBGroup
        if ismember([n1, n2], NB{j}, 'rows') || ismember([n2, n1], NB{j}, 'rows')
            bgroup = j;
            break
        end
    end
    B2E(i, :) = [elem, face, bgroup];
    %find the vector parallel to the boundary and counterclockwise to elem, 
    %then rotate it by 90 degrees clockwise
    idx1 = mod(face, 3) + 1; idx2 = mod(face + 1, 3) + 1;
    faceVec = p(:, E2N(elem, idx2)) - p(:, E2N(elem, idx1));
    R = [0, 1; -1, 0];        %rotation matrix
    nVec = R*faceVec;
    nVec = nVec/norm(nVec);   %normal vector pointing outward from elem
    Bn(i, :) = nVec';
end

%Area matrix
Area = zeros(nElemTot, 1);
for i = 1 : nElemTot
    l1 = p(:, E2N(i, 2)) - p(:, E2N(i, 1));
    l2 = p(:, E2N(i, 3)) - p(:, E2N(i, 1));
    Area(i) = 0.5*abs(l1(1)*l2(2) - l1(2)*l2(1));  %magnitude of the cross product
end

end