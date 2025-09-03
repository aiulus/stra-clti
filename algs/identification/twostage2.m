function out = twostage2(Y, Ts, nord, rough_pen)
%TWOSTAGE2 Functional two-stage estimator using smoothing splines.
%   Attempts to mimic the R/fda pipeline; auto-fallbacks described below.
%
%   Backends (auto-detected):
%     1) fdaM (create_bspline_basis/fdPar/smooth_basis/inprod/deriv)
%        -> Builds y2cMap, Jb, Dmat to produce S and L exactly as in R.
%     2) csaps-based operator construction:
%        -> Builds linear operators H (smoothing) and D (smoothed deriv.)
%           so that S = H'*W*H and L = D'*W*H with trapezoidal W.
%
%   Inputs:
%     Y        : d x n
%     Ts       : 1 x n (in ascending order)
%     nord     : spline order (default 4 == cubic)
%     rough_pen: roughness penalty (default 1e-3)
%
%   Outputs:
%     Ahat, S, L, xt_hat (struct with fields t,yhat,pp per state), x0

if nargin < 3, nord = 4; end
if nargin < 4, rough_pen = 1e-3; end
[d,n] = size(Y);
Ts = Ts(:).';
if numel(Ts) ~= n, error('Length(Ts) must equal size(Y,2).'); end

% Quadrature weights (trapezoidal)
w = trapz_weights(Ts);               % 1 x n
W = diag(w);

% Try FDA-M style pipeline
hasFDA = (exist('create_bspline_basis','file')==2) && ...
         (exist('fdPar','file')==2) && ...
         (exist('smooth_basis','file')==2) && ...
         (exist('inprod','file')==2) && ...
         (exist('fd','file')==2) && ...
         (exist('deriv','file')==2);

if hasFDA
    rangeval = [Ts(1) Ts(end)];
    nbasis   = n + nord - 2;
    basis    = create_bspline_basis(rangeval, nbasis, nord, Ts);
    % build y2cMap by smoothing the identity columns
    y2cMap = zeros(n, nbasis);
    harmaccelLfd = int2Lfd(2);
    fdParobj = fdPar(basis, harmaccelLfd, rough_pen);
    I = eye(n);
    for j = 1:n
        fdobj = smooth_basis(Ts, I(:,j).', fdParobj); % fit basis to delta column
        y2cMap(j,:) = getcoef(fdobj).';               % coefficients for column j
    end
    % inner products of basis
    Jb = inprod(basis, basis);
    % derivative operator matrix in coefficient space
    Dfd = deriv(fd(eye(nbasis), basis));              % derivative of basis coefficients
    Dmat = getcoef(Dfd);                              % nbasis x nbasis (maps coeffs of f to coeffs of f')

    S = y2cMap.' * Jb * y2cMap;
    L = y2cMap.' * Dmat.' * Jb * y2cMap;

    % smoothed curves for return
    xt_hat = fsmooth(Y, Ts, rough_pen, nord, false, 5);
else
    % csaps-based operator construction (linear operators)
    hasCSAPS = exist('csaps','file')==2 && exist('fnder','file')==2 && exist('fnval','file')==2;
    p = 1 / (1 + max(rough_pen, 0));
    % Build H (n x n): maps raw samples -> smoothed values at Ts
    % Build D (n x n): maps raw samples -> smoothed derivative at Ts
    H = zeros(n,n); D = zeros(n,n);
    E = eye(n);
    for j = 1:n
        yj = E(:,j);
        if hasCSAPS
            pp = csaps(Ts, yj.', p);
            H(:,j) = fnval(pp, Ts).';
            ppd    = fnder(pp, 1);
            D(:,j) = fnval(ppd, Ts).';
        else
            % Fallback: same as fsmooth fallback + finite-diff on smoothed
            % (keeps linearity)
            alpha = max(rough_pen, 1e-3);
            D2 = diff(eye(n),2);
            A = speye(n) + alpha*(D2.'*D2);
            z = A \ yj;          % smoothed samples
            H(:,j) = z;
            % derivative via centered difference on smoothed data
            Dz = zeros(n,1);
            dt = diff(Ts);
            % interior points
            for k = 2:n-1
                Dz(k) = (z(k+1)-z(k-1)) / (Ts(k+1)-Ts(k-1));
            end
            % endpoints (forward/backward)
            Dz(1) = (z(2)-z(1)) / (Ts(2)-Ts(1));
            Dz(n) = (z(n)-z(n-1)) / (Ts(n)-Ts(n-1));
            D(:,j) = Dz;
        end
    end

    S = H.' * W * H;      % ~~ \int xhat xhat^T
    L = D.' * W * H;      % ~~ \int xhat' xhat^T

    % smoothed curves for return (reuse fsmooth)
    xt_hat = fsmooth(Y, Ts, rough_pen, nord, false, 5);
end

Ahat = Y * L * Y.' * rsolve(Y * S * Y.');
out  = struct('Ahat',Ahat,'S',S,'L',L,'xt_hat',xt_hat,'x0',Y(:,1).');

end