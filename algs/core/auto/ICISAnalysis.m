function out = ICISAnalysis(A, x0, n_digits)
%ICISANALYSIS Identifiability at initial condition x0 (ICIS & related).
%   OUT = ICISANALYSIS(A, x0, N_DIGITS)
%
%   OUT fields:
%     Identifiable, Ident1, Ident2, ICIS, w0k_norm, I0, Iplus,
%     Jordan (struct: J,Qmat,Qinv,K1,K2,Lgap,RepEigenIdx)

%% TODO: Add references

if nargin < 3, n_digits = 3; end
[p1,p2] = size(A);
if p1 ~= p2, error('A must be square'); end
if numel(x0) ~= p1, error('Length of x0 must equal size(A,1)'); end

% 1) Jordan over R
jcf   = realJordan(A);
Qinv  = jcf.Qinv;
K1    = jcf.K1;  K2 = jcf.K2;  K = K1 + K2;

% 2) Repeated eigenvalues (rounded to n_digits)
lmb  = eig(A);
fac  = 10^n_digits;
lmbR = round(real(lmb)*fac)/fac + 1i*round(imag(lmb)*fac)/fac;
% detect repeats after sorting for stability
[ls,ord] = sortrows([real(lmbR), imag(lmbR)]);
runs = 1;
RepEigenIdx = [];
if ~isempty(ls)
    cnt = 1;
    blocks = {};
    for i = 2:size(ls,1)
        if all(abs(ls(i,:) - ls(i-1,:)) < eps)
            cnt = cnt + 1;
        else
            blocks{end+1} = cnt; %#ok<AGROW>
            cnt = 1;
        end
    end
    blocks{end+1} = cnt;
    runs = cell2mat(blocks);
end
if max(runs) == 1
    Ident1 = true;
else
    Ident1 = false;
    % Return first repeated block's original indices (best-effort)
    r0 = find(runs > 1, 1, 'first');
    s  = sum(runs(1:r0-1))+1; e = s + runs(r0) - 1;
    RepEigenIdx = ord(s:e);
end
Lgap = Lgapfun(lmb);

% 3) |w0,k|
xtilde0 = Qinv * x0(:);
w0k_norm = zeros(K,1);
if K1 > 0
    w0k_norm(1:K1) = round(abs(xtilde0(1:K1)), n_digits);
end
if K2 > 0
    for k = 1:K2
        k2 = K1 + 2*k - 1;
        w0k_norm(K1+k) = round( sqrt(xtilde0(k2)^2 + xtilde0(k2+1)^2), n_digits);
    end
end
ICIS   = min(w0k_norm);
Ident2 = (ICIS > 0);
ident  = Ident1 && Ident2;

% 5) I0 and Iplus masks lifted to full state dimension
a0 = double(w0k_norm <= 0); % 1 if unexcited subspace
if K1 > 0, aa_real = a0(1:K1);
else,      aa_real = []; end
if K2 > 0, aa_comp = repelem(a0(K1+1:K), 2).';
else,      aa_comp = []; end
aa    = [aa_real(:); aa_comp(:)];
I0    = diag(aa);
Iplus = eye(p1) - I0;

out = struct( ...
    'Identifiable', ident, ...
    'Ident1', Ident1, ...
    'Ident2', Ident2, ...
    'ICIS', ICIS, ...
    'w0k_norm', w0k_norm, ...
    'I0', I0, ...
    'Iplus', Iplus, ...
    'Jordan', struct('J',jcf.J,'Qmat',jcf.Qmat,'Qinv',jcf.Qinv, ...
                     'K1',K1,'K2',K2,'Lgap',Lgap,'RepEigenIdx',RepEigenIdx) );
end
