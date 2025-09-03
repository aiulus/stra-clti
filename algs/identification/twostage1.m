function out = twostage1(Y, tstep)
%TWOSTAGE1 Two-stage estimator using finite differences (simple).
%   out.Ahat = Y*L*Y' / (Y*Y'), with L from first-order difference / tstep.
%
%   Inputs:
%     Y     : d x n matrix (rows = states, cols = timepoints)
%     tstep : constant time step (scalar)
%
%   Outputs:
%     Ahat, S (n x n identity), L (n x n difference), x0 (1 x d)

[d,n] = size(Y);
if n < 2, error('Need at least 2 time points.'); end
if ~isscalar(tstep) || tstep <= 0, error('tstep must be positive scalar'); end

% L = [D  0] where D is forward difference along time / tstep
D = -[eye(n-1); zeros(1,n-1)] + [zeros(1,n-1); eye(n-1)];
L = [D, zeros(n,1)] / tstep;

Ahat = Y * L * Y.' * rsolve(Y*Y.');
out = struct('Ahat',Ahat, 'S',eye(n), 'L',L, 'x0', Y(:,1).');
end
