%Basis functions evaluated at given points in reference space
function phi = basis(xiEta, p, coeff)

n = size(xiEta, 2);
np = (p + 1)*(p + 2)/2;
funs = zeros(n, np);

for i = 1 : n
    k = 1;
    xi = xiEta(1, i); eta = xiEta(2, i);
    for s = 0 : p
        for r = 0 : p - s
            funs(i, k) = xi^r*eta^s;
            k = k + 1;
        end
    end
end

phi = funs*coeff;  %(n x np)*(np x np) matrix, each entry is one basis 
                   %function evaluated at one point

end