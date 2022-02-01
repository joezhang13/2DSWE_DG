function mesh = readgri(grifile)
% function mesh = readgri(grifile);
%
% This function reads a text .gri file into a data structure called
% mesh.
%
% INPUTS: 
%   grifile = name of .gri file
%
% OUTPUTS:
%   mesh = data structure:
%    mesh.Dim   = dimension
%    mesh.nNode = number of nodes
%    mesh.Node  = nNode x dim array of node coordinates
%    mesh.nElem  = number of elements
%    mesh.B = boundary information
%         B.nbfgrp = number of boundary face groups
%         B.nbface = [nbfgrp] number of faces in each group
%         B.nnode = [nbfgrp] number of nodes per face in each group
%         B.title = [nbfgrp] title of each group
%         B.nodes = cell array of boundary nodes
%    mesh.QBasis = basis type for elements
%    mesh.QOrder = order type for elements
%    mesh.Elem   = array of element node numbers (1-based)
%

% initialize mesh structure
mesh = struct;

% open the file for writing
fid = fopen(grifile, 'r');

% Read in nodes
A = fscanf(fid,'%d', 3);
mesh.nNode    = A(1);
nelemtot = A(2);
mesh.Dim = A(3);
mesh.Node = zeros(mesh.nNode, mesh.Dim);
for inode = 1:mesh.nNode
  A = fscanf(fid, '%lf', mesh.Dim);
  mesh.Node(inode,:) = A(1:mesh.Dim)';
end

% Read boundary info
A = fscanf(fid, '%d', 1);
B.nbfgrp = A(1);
B.nbface = zeros(B.nbfgrp,1);
B.nnode  = zeros(B.nbfgrp,1);
B.title  = cell(B.nbfgrp,1);
B.nodes  = cell(B.nbfgrp,1);
for ibfgrp = 1:B.nbfgrp
  fgets(fid);
  sline = fgets(fid);
  [B.nbface(ibfgrp), B.nnode(ibfgrp), B.title(ibfgrp)] = strread(sline, '%d %d %s');
  N = zeros(B.nbface(ibfgrp), B.nnode(ibfgrp));
  for ibface = 1:B.nbface(ibfgrp)
    A = fscanf(fid, '%d', B.nnode(ibfgrp));
    N(ibface,:) = A';
  end
  B.nodes{ibfgrp} = N;
end
mesh.B = B;

% Read in elements 
fgets(fid);
sline = fgets(fid);
[nelem, p, sbasis] = strread(sline, '%d %d %s');
mesh.nElem = nelem;
mesh.QBasis = sbasis{1};
mesh.QOrder = p;
switch sbasis{1}
  case 'TriLagrange'
    nnode = (p+1)*(p+2)/2;
  otherwise
    error('element type not understood');
end

E = zeros(nelem, nnode);
  
for elem = 1:nelem
  E(elem,:) = fscanf(fid, '%d', nnode);  % HO nodes included
end

mesh.Elem = E;

% close file
fclose(fid);
