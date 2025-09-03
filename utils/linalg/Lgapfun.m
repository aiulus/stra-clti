function Lgap = Lgapfun(lambdas)
%LGAPFUN Minimum squared distance between (complex) eigenvalues.
%   Lgap = min_{i<j} |lambda_i - lambda_j|^2

lambdas = lambdas(:);
z = lambdas;
% pairwise distances without forming giant matrix
Lgap = inf;
for i = 1:numel(z)-1
    dz = z(i+1:end) - z(i);
    d2 = min(abs(dz).^2);
    if d2 < Lgap, Lgap = d2; end
end
if isinf(Lgap), Lgap = 0; end
end
