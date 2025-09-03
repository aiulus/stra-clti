function X = rsolve(mat, d_prop, dmin, dmax)
%RSOLVE Numerically robust inverse using SVD thresholding.
%   X = RSOLVE(MAT, d_prop=1e-6, dmin=1e-9, dmax=1e9)

if nargin < 2, d_prop = 1e-6; end
if nargin < 3, dmin   = 1e-9;  end
if nargin < 4, dmax   = 1e9;   end

[U,S,V] = svd(mat, 'econ');
s = diag(S);
thresh = min(max(sum(s) * d_prop, dmin), dmax);
sInv = 1 ./ max(s, thresh);
X = V * diag(sInv) * U.';
end
