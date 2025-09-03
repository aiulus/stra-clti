function out = realJordan(A, tol)
%REALJORDAN Jordan-like real block form and real eigenbasis (numeric).
%   OUT = JORDANREAL(A, TOL) computes a real block-diagonal J and a real
%   basis Q such that A ≈ Q*J/Q (no nilpotent blocks supported).
%
%   Returned struct OUT has fields:
%     J        - block-diagonal: 1x1 real lambdas; 2x2 [a b; -b a] blocks
%     Qmat     - real "generalized" eigenvectors (columns)
%     Qinv     - numerically robust inverse of Qmat via SVD thresholding
%     lambdas  - eigenvalues of A (complex, as returned by eig)
%     K1       - # of real eigenvalues (|Imag| < tol)
%     K2       - # of complex-conjugate pairs (Imag > tol counts)
%
%   Notes:
%   * Not numerically stable for difficult matrices; assumes no nilpotent
%     cells in the true JCF.
%   * Pairs complex conjugates by selecting the positive-imag partner.
%
%   See also: eig, trace, rsolve

if nargin < 2, tol = 1e-6; end
[n,m] = size(A);
if n ~= m, error('A must be square.'); end

[V,D] = eig(A);
lambdas = diag(D);
imPart  = imag(lambdas);

isReal  = abs(imPart) < tol;
realIdx = find(isReal);
posIdx  = find(imPart >  tol);  % only the positive-imag partners

K1 = numel(realIdx);
K2 = numel(posIdx);

J    = zeros(n);
Qmat = zeros(n);

% Place real eigenvalues/vecs first
for k = 1:K1
    idx = realIdx(k);
    J(k,k)     = real(lambdas(idx));
    Qmat(:,k)  = real(V(:,idx));
end

% Complex pairs: use the +imag partner's eigenvector v
for k = 1:K2
    idx = posIdx(k);
    a   = real(lambdas(idx));
    b   = imag(lambdas(idx)); % > 0 by construction
    k2  = K1 + 2*k - 1;

    J(k2:k2+1, k2:k2+1) = [a  b; -b  a];
    v  = V(:,idx);
    Qmat(:,k2)   = real(v);
    Qmat(:,k2+1) = imag(v);
end

Qinv = rsolve(Qmat);
out  = struct('J',J,'Qmat',Qmat,'Qinv',Qinv, ...
              'lambdas',lambdas,'K1',K1,'K2',K2);
end
