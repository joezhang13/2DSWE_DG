%Initial condition
function [U0, dt] = initial(mesh, CFL0, p)

node = mesh.Node;
E2N = mesh.Elem;
ne = mesh.nElem;
np = (p + 1)*(p + 2)/2;
[~, pref, ~] = TriLagrange2D(p);  %Lagrange nodes in the reference space

U0 = zeros(ne*np, 3);
pn = zeros(ne*np, 2);          %coordinates of the Lagrange nodes
d = zeros(ne, 1);              %perimeter of each cell
for i = 1 : ne
    pp = node(E2N(i, :), :);
    l1 = pp(3, :) - pp(2, :);
    l2 = pp(1, :) - pp(3, :);
    l3 = pp(2, :) - pp(1, :);
    A = 0.5*abs(l1(1)*l2(2) - l2(1)*l1(2));      %cell area
    P = norm(l1) + norm(l2) + norm(l3);          %cell perimeter
    d(i) = 2*A/P;                                %measure of the cell size
    x1 = pp(1, :)'; x2 = pp(2, :)'; x3 = pp(3, :)';
    J = [x2(1)-x1(1), x3(1)-x1(1); x2(2)-x1(2), x3(2)-x1(2)];
    if p == 0
        xn = (x1 + x2 + x3)/3;
    else
        xn = x1 + J*pref;                            %Lagrange nodes in global space
    end
    pn((i-1)*np+1 : i*np, :) = xn';
end
x = pn(:, 1); y = pn(:, 2);
h = 1 + 0.3*exp(-50*(x - 1.5).^2 - 50*(y - 0.7).^2);
U0(:, 1) = h;

g = 9.8;
smax = max((g*h).^0.5);
dmin = min(d);
dt = CFL0/(p + 1)*dmin/smax;

end