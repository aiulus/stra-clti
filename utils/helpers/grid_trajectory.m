function X = grid_trajectory(A, x0, tgrid)
%TRAJECTORY Compute x(t) = expm(A t) x0 for each t in tgrid.
[V,D] = eig(A);
Vinvt = inv(V);
lam   = diag(D);
X     = zeros(numel(tgrid), numel(x0));
for i = 1:numel(tgrid)
    ti = tgrid(i);
    M  = V * diag(exp(lam*ti)) * Vinvt;
    xi = real(M * x0(:));  % discard tiny imag from numerics
    X(i,:) = xi(:).';
end
end