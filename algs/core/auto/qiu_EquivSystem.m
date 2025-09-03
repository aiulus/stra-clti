function Atilde = qiu_EquivSystem(A, x0, Dmat, varargin)
%EQUIVSYSTEM Build an equivalent system matrix for (A, x0).
%   Atilde = A + Q*I0*D*I0*Q^{-1} where I0 encodes unexcited subspaces.
%
%   When A is identifiable at x0, Atilde == A for all Dmat. When Dmat==0,
%   Atilde == A.

if ~isequal(size(A), size(Dmat))
    error('A and Dmat must have equal dimensions.');
end
[p1,p2] = size(A);
if p1 ~= p2, error('A must be square'); end
if numel(x0) ~= p1, error('Length of x0 must equal size(A,1)'); end

obj    = ICISAnalysis(A, x0, varargin{:});
Qmat = obj.Jordan.Qmat;
Qinv = obj.Jordan.Qinv;
I0   = obj.I0;

Atilde = A + Qmat * (I0 * (Dmat * (I0 * Qinv)));
end
