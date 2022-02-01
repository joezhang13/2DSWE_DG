function [IE,BE] = edgehash(E2N)
% This function identifies interior and boundary edges, and their
% connectivities, in a triangular mesh given an element-to-node array.
%
% INPUT : E2N = [nelem x 3] array mapping triangles to nodes
% OUTPUT: IE = [niedge x 4] array giving (n1, n2, elem1, elem2)
%              information for each interior edge
%         BE = [nbedge x 3] array giving (n1, n2, elem)
%              information for each boundary edge

nelem = size(E2N,1);            % number of elements
nnode = max(max(E2N));          % number of nodes
H = sparse(nnode,nnode);        % Create a hash list to identify edges
IE = zeros(ceil(nelem*3/2), 4); % (over) allocate interior edge array
niedge = 0;                     % number of interior edges (running total)

% Loop over elements and identify all edges
for elem = 1:nelem
  nv = E2N(elem,1:3);
  for edge = 1:3
    n1 = nv(mod(edge  ,3)+1);
    n2 = nv(mod(edge+1,3)+1);
    if (H(n1,n2) == 0) % edge hit for the first time
      % could be a boundary or interior; assume boundary
      H(n1,n2) = elem;  H(n2,n1) = elem;
    else % this is an interior edge, hit for the second time
      oldelem = H(n1,n2);
      if (oldelem < 0), error 'Mesh input error'; end;
      niedge = niedge+1;
      IE(niedge,:) = [n1,n2, oldelem, elem];
      H(n1,n2) = -1;  H(n2,n1) = -1;
    end
  end
end

IE = IE(1:niedge,:);  % clip IE

% find boundary edges
[I,J] = find(triu(H)>0);  
BE = [I, J, zeros(size(I))];
for b = 1:size(I,1), BE(b, 3) = H(I(b),J(b)); end;


