%Solution plotting
function solplot(mesh, U, p, t)

%mesh information and solution order
node = mesh.Node;
E2N = mesh.Elem;
ne = mesh.nElem;
np = (p + 1)*(p + 2)/2;

%pre-computed plotting points on the reference triangle
[~, pplot, ~] = TriLagrange2D(p + 1);     %coordinates of plotting points
[coeff, ~, ~] = TriLagrange2D(p);
phi = basis(pplot, p, coeff);             %interpolation matrix (from Lagrange nodes to plotting points)
nsub = (p + 1)^2;                         %number of sub-elements
subE2N = zeros(nsub, 3);                  %element-to-node matrix for the sub-elements
elem = 1; pidx = 1;                       %store the sub-element index and the point index
for i = 1 : p + 1
    subE2N(elem, :) = [pidx, pidx+1, pidx+p+3-i];
    elem = elem + 1; pidx = pidx + 1;
    for j = 2 : p + 2 - i
        subE2N(elem, :) = [pidx, pidx+p+3-i, pidx+p+2-i];
        subE2N(elem+1, :) = [pidx, pidx+1, pidx+p+3-i];
        elem = elem + 2; pidx = pidx + 1;
    end
    pidx = pidx + 1;
end

%Loop over elements for plotting
figure;
for i = 1 : ne
    indice = (i-1)*np+1 : i*np;
    h = U(indice, 1);
    h = phi*h;
    pp = node(E2N(i, :), :);
    x1 = pp(1, :)'; x2 = pp(2, :)'; x3 = pp(3, :)';
    J = [x2(1)-x1(1), x3(1)-x1(1); x2(2)-x1(2), x3(2)-x1(2)];          %mapping Jacobian
    x = x1 + J*pplot; x = x';
    patch('Faces',subE2N,'Vertices',x,'FaceVertexCData',h,'FaceColor','interp','EdgeColor','none');
    hold on
end
axis equal
xlim([0 2]);
ylim([0 1]);
colorbar;
caxis([0.9 1.1]);
if t == 0
    caxis([0.7 1.3]);
end
set(gca, 'FontSize', 15);

end