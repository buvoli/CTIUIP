function [C, FD] = iPolyCoefficients(z, b, numeric_type)
% IPOLYCOEFFICIENTS: computes the coefficients for integrating a polynomial 
% p(x) from x=a to x=b. The polynomial p(x) interpolates data  
% values specified  points. 
%
% Parameters
%   z   - (vector) node locations
%   b   - (matrix) b(:,1) are integration start-pint and b(:,2) are integration endpoints
% Return Values
%   C   - matrix of weights for computing outputs
%   FD  - associated finite difference matrix for computing derivatives of
%         p(x) at x = 0.
% Example: iPolyCoefficients([-1 0 1], [0 1])

q = length(z);
num_ints = size(b,1);
if(nargin == 2)
    numeric_type = 'vpa';
end

if(strcmp(numeric_type, 'sym'))
    I = diag(sym(1)./sym(1:q));
    M = sym(zeros(q));
    E = sym(zeros(num_ints, q));
elseif(strcmp(numeric_type, 'vpa'))
    I = diag(vpa(1)./vpa(1:q));
    M = vpa(zeros(q));
    E = vpa(zeros(num_ints, q));
else
    I = diag(1./(1:q));
    M = zeros(q);
    E = zeros(num_ints, q);
end

for i=1:q
    for j=1:q
        M(i,j) = z(i)^(j-1);
    end
end


% -- Evaluation --------------
if(nargin >=  2)
    for i=1:num_ints
        for j=1:q
            E(i,j) = b(i,2)^j - b(i,1)^j;
        end
    end
    C = transpose(transpose(M) \ transpose(E*I)); % equivalent and more stable than C = E*I*inv(M);
else
    C = inv(M);
end
if(nargout == 2)
    FD = inv(M);
end
end