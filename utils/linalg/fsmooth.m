function Xt = fsmooth(Y, Ts, rough_pen, norder, doplot, n_plot)
%FSMOOTH Simple wrapper for smoothing splines (multivariate).
%   Xt is a struct array of smoothing spline fits (per row of Y).
%
%   Priority of backends:
%     1) fdaM (create_bspline_basis/fdPar/smooth_basis) if available
%     2) csaps (Curve Fitting Toolbox)
%     3) Fallback: cubic smoothing via piecewise cubic Hermite + mild
%        Tikhonov regularization (less accurate)

if nargin < 3, rough_pen = 1e-3; end
if nargin < 4, norder    = 4;    end
if nargin < 5, doplot    = false; end
if nargin < 6, n_plot    = 5;    end

[d,n] = size(Y);
Ts = Ts(:).';
if numel(Ts) ~= n, error('Length(Ts) must equal size(Y,2).'); end

Xt = struct('t',Ts,'yhat',[],'pp',[]); Xt = repmat(Xt, d, 1);

% --- Backend 1: fdaM if present ---
hasFDA = (exist('create_bspline_basis','file')==2) && ...
         (exist('fdPar','file')==2) && (exist('smooth_basis','file')==2);

if hasFDA
    rangeval = [Ts(1), Ts(end)];
    nbasis   = n + norder - 2;
    basisobj = create_bspline_basis(rangeval, nbasis, norder, Ts);
    % lambda ~ rough_pen (user provided)
    harmaccelLfd = int2Lfd(2);
    fdParobj     = fdPar(basisobj, harmaccelLfd, rough_pen);
    for i = 1:d
        fdobj = smooth_basis(Ts, Y(i,:), fdParobj);
        Xt(i).pp   = fdobj; % fdaM fd object
        Xt(i).yhat = eval_fd(Ts, fdobj);
    end
else
    % --- Backend 2: csaps if available; else Fallback ---
    hasCSAPS = exist('csaps','file')==2;
    % Map rough_pen (lambda-like) to csaps p \in [0,1]: p ~ 1/(1+lambda)
    p = 1 / (1 + max(rough_pen, 0));

    for i = 1:d
        yi = Y(i,:);
        if hasCSAPS
            pp = csaps(Ts, yi, p);
            yhat = fnval(pp, Ts);
            Xt(i).pp   = pp;
            Xt(i).yhat = yhat;
        else
            % Fallback: mild smoothing via ridge on finite-diff operator
            % Solve min ||y - z||^2 + alpha ||D z||^2
            alpha = max(rough_pen, 1e-3);
            D = diff(eye(n),2);        % second-difference
            A = speye(n) + alpha*(D.'*D);
            z = A \ yi(:);
            yhat = z(:).';
            pp = spline(Ts, yhat);     % store as a piecewise cubic
            Xt(i).pp   = pp;
            Xt(i).yhat = yhat;
        end
    end
end

if doplot
    n2 = min(d, n_plot);
    figure; hold on;
    plot(Ts, Y(1:n2,:).', '.', 'DisplayName', 'data');
    for i = 1:n2
        plot(Ts, Xt(i).yhat, 'LineWidth', 1.25, 'DisplayName', sprintf('fit %d',i));
    end
    xlabel('t'); ylabel('x(t)'); title('Smoothing splines');
    grid on; hold off;
end
end
