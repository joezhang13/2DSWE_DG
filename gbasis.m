%The gradients of the basis functions evaluated at given points in reference space
function dphidxi = gbasis(xiEta, p, coeff)

n = size(xiEta, 2);
np = (p + 1)*(p + 2)/2;
gfuns = zeros(n, np, 2);

for i = 1 : n
    k = 1;
    xi = xiEta(1, i); eta = xiEta(2, i);
    for s = 0 : p
        for r = 0 : p - s
            if r == 0
                gfuns(i, k, 1) = 0;
            else
                gfuns(i, k, 1) = r*xi^(r - 1)*eta^s;
            end
            if s == 0
                gfuns(i, k, 2) = 0;
            else
                gfuns(i, k, 2) = xi^r*s*eta^(s - 1);
            end
            k = k + 1;
        end
    end
end

dphidxi(:, :, 1) = gfuns(:, :, 1)*coeff;
dphidxi(:, :, 2) = gfuns(:, :, 2)*coeff;

end