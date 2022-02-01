%Calculate coefficients for full-order Lagrange basis of order p and the
%coordinates of the Lagrange nodes in the reference space (reference
%element is a unit isoceles right triangle).
function [C, refNodes, faceIdx] = TriLagrange2D(p)

xi  = linspace(0,1,p+1);  eta = xi;
N = (p+1)*(p+2)/2; % number of basis functions
A = zeros(N,N);
i = 1; % build A-matrix
for iy=0:p
    for ix=0:p-iy  % loop over nodes
        k = 1;
        for s=0:p
            for r=0:p-s % loop over monomials
                A(i,k) = xi(ix+1)^r * eta(iy+1)^s;
                k = k+1;
            end
        end
        i = i + 1;
    end
end
C = inv(A);

%obtain the coordinates of the Lagrange nodes in the reference space
refNodes = zeros(2, N);
k = 1;
for i = 0 : p 
    for j = 0 : p - i
        xiref = xi(j + 1);
        etaref = eta(i + 1);
        refNodes(:, k) = [xiref; etaref];
        k = k + 1;
    end
end

%obtain the indice of the Lagrange nodes on the cell faces in the
%reference space (counterclockwise)
faceIdx = zeros(3, p + 1);
i = 1 : p + 1;
faceIdx(1, :) = i.*(2*p+3-i)/2;
faceIdx(2, :) = (i-1).*(2*p+4-i)/2+1;
faceIdx(3, :) = i;

end

% % print out coefficients 
% for s=0:p
%     for r=0:p-s 
%         fprintf('%11s', sprintf('x^%d*y^%d', r, s));
%     end
% end
% fprintf('\n'); fprintf(strcat(repmat(' %10g', 1,N), '\n'), C);
