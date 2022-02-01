%Roe flux for 2D Shallow water equation
function [F, smax] = RoeFlux(UL,UR,n)

g = 9.8;
nn = size(UL, 2);

%left state
hL = UL(1, :);
uL = UL(2, :)./hL;
vL = UL(3, :)./hL;
for j = 1 : nn
    if hL(j) < 0
        disp('Non-physical L state!');
        disp(UL(:, j));
    end
end
unL = uL.*n(1, :) + vL.*n(2, :);
%left flux
FL = zeros(3, nn);
FL(1, :) = hL.*unL;
FL(2, :) = hL.*uL.*unL + 0.5*g*hL.^2.*n(1, :);
FL(3, :) = hL.*vL.*unL + 0.5*g*hL.^2.*n(2, :);

%right state
hR = UR(1, :);
uR = UR(2, :)./hR;
vR = UR(3, :)./hR;
for j = 1 : nn
    if hR(j) < 0
        disp('Non-physical R state!');
        disp(UR(:, j));
    end
end
unR = uR.*n(1, :) + vR.*n(2, :);
%right flux
FR = zeros(3, nn);
FR(1, :) = hR.*unR;
FR(2, :) = hR.*uR.*unR + 0.5*g*hR.^2.*n(1, :);
FR(3, :) = hR.*vR.*unR + 0.5*g*hR.^2.*n(2, :);

%difference in states
dh = hR - hL;
dhu = UR(2, :) - UL(2, :);
dhv = UR(3, :) - UL(3, :);

%Roe average (arithmetic average for SWE)
ha = 0.5*(hL + hR);
ua = 0.5*(uL + uR);
va = 0.5*(vL + vR);
c = (g*ha).^0.5;
un = ua.*n(1, :) + va.*n(2, :);

%eigen values
lambda = zeros(3, nn);
lambda(1, :) = un;
lambda(2, :) = un - c;
lambda(3, :) = un + c;

%entropy fix
epsilon = 0.01*c;
for j = 1 : nn
    for i = 1 : 3
        if abs(lambda(i, j)) < epsilon(j)
            lambda(i, j) = 0.5*(epsilon(j) + lambda(i, j)^2/epsilon(j));
        end
    end
end

lambda = abs(lambda);

%coefficients in the stabilization term
s1 = 0.5*(lambda(2, :) + lambda(3, :));
s2 = 0.5*(lambda(2, :) - lambda(3, :));
A1 = s2./c.*un + s1;
B1 = -s2./c.*(dhu.*n(1, :) + dhv.*n(2, :));
A2 = lambda(1, :);
B2 = (A1 - A2).*dh + B1;
C2 = ((lambda(1, :) - s1).*un - s2.*c).*dh + (s1 - lambda(1, :)).*(dhu.*n(1, :) + dhv.*n(2, :));

%flux assembly
F = zeros(3, nn);
F(1, :) = 0.5*(FL(1, :) + FR(1, :)) - 0.5*(A1.*dh  + B1);
F(2, :) = 0.5*(FL(2, :) + FR(2, :)) - 0.5*(A2.*dhu + B2.*ua + C2.*n(1, :));
F(3, :) = 0.5*(FL(3, :) + FR(3, :)) - 0.5*(A2.*dhv + B2.*va + C2.*n(2, :));

smax = max(lambda, [], 1);

end