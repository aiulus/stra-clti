function[figA, figAt] = example_Ax0_identifiability(varargin)
% PLOT_FIVE_LAYOUT_MODAL_BASIS  Five-plot (1,2,2) layouts for A and Ã.
%
%   Builds Q from a real modal/Jordan-like basis: one real eigenvector (L1)
%   and the real/imag parts of a complex eigenvector (L2 plane). Shows:
%     Row 1: subspaces only (L1 line, L2 plane)
%     Row 2: projections of x0^(a) [left] and x0^(b) [right] onto L1/L2
%     Row 3: trajectories from x0^(a) [left] and x0^(b) [right]
%
%   By default uses the problem data embedded below and saves:
%     figs/A_layout.png
%     figs/Atilde_layout.png
%
%   [figA, figAt] = plot_five_layout_modal_basis()                 % default
%   [..] = plot_five_layout_modal_basis('A',A,'Atilde',At,...)     % override
%
%   Name-Value options:
%     'A'            : 3x3 real matrix (default below)
%     'Atilde'       : 3x3 real matrix (default below)
%     'z0_a'         : 3x1 modal coords for A (default [2;-1;0])
%     'z0_b'         : 3x1 modal coords for A (default [0;-2;3])
%     'T'            : horizon (default 6.0)
%     'N'            : time samples (default 600)
%     'OutDir'       : output directory (default 'figs')
%     'SaveA'        : filename for A figure (default 'A_layout.png')
%     'SaveAt'       : filename for Ã figure (default 'Atilde_layout.png')
%     'Show'         : logical, show figures (default true)
%     'Span'         : subspace drawing span (default 2.4)
%     'PlaneAlpha'   : plane transparency (default 0.22)
%     'UseEigExpm'   : use eig-based expm (default false; uses expm otherwise)
%
%   Returns figure handles figA and figAt.

% -------- Defaults (problem data) -----------------------------------------
A_def = [ 0  1 -1;
          2  0  0;
          3  1  0];

At_def = [ 1 -1  0;
           0  4 -2;
           2  3 -1];

z0a_def = [ 2; -1; 0];
z0b_def = [ 0; -2; 3];

% -------- Parse inputs ----------------------------------------------------
p = inputParser;
p.addParameter('A',        A_def,      @(x)isreal(x)&&isequal(size(x),[3,3]));
p.addParameter('Atilde',   At_def,     @(x)isreal(x)&&isequal(size(x),[3,3]));
p.addParameter('z0_a',     z0a_def,    @(x)isreal(x)&&isequal(size(x),[3,1]));
p.addParameter('z0_b',     z0b_def,    @(x)isreal(x)&&isequal(size(x),[3,1]));
p.addParameter('T',        6.0,        @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('N',        600,        @(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('OutDir',   'figs',     @(s)ischar(s)||isstring(s));
p.addParameter('SaveA',    'A_layout.png',      @(s)ischar(s)||isstring(s));
p.addParameter('SaveAt',   'Atilde_layout.png', @(s)ischar(s)||isstring(s));
p.addParameter('Show',     true,       @(x)islogical(x)&&isscalar(x));
p.addParameter('Span',     2.4,        @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('PlaneAlpha',0.22,      @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('UseEigExpm',false,     @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
S = p.Results;

% -------- Output directory ------------------------------------------------
if ~exist(S.OutDir,'dir'), mkdir(S.OutDir); end

% -------- Build real modal bases -----------------------------------------
[Q_A,   q1_A,   U2_A]  = real_modal_basis(S.A);
[Q_At,  q1_At,  U2_At] = real_modal_basis(S.Atilde);

% Initial conditions in state coords defined via A's modal basis
x0_a = Q_A * S.z0_a;
x0_b = Q_A * S.z0_b;

% -------- Build figures ---------------------------------------------------
tagA   = 'A';
tagAt  = 'Ã';   % display tag only

figA  = make_figure_for_system(S.A,    Q_A,  q1_A,  U2_A,  x0_a, x0_b, tagA,  ...
                               S.T, S.N, fullfile(S.OutDir,S.SaveA), S.Show, ...
                               S.Span, S.PlaneAlpha, S.UseEigExpm);

figAt = make_figure_for_system(S.Atilde,Q_At, q1_At, U2_At, x0_a, x0_b, tagAt, ...
                               S.T, S.N, fullfile(S.OutDir,S.SaveAt), S.Show, ...
                               S.Span, S.PlaneAlpha, S.UseEigExpm);
end

% ========================= Helper functions ===============================

function [Q, q1, U2] = real_modal_basis(A, tol)
% Build a real modal basis Q for a 3x3 real A with one real eigenvalue and
% a complex-conjugate pair. q1 is the real eigenvector (L1), U2 = [u2 v2]
% spans the invariant plane (L2) using Re/Im of a complex eigenvector.
if nargin < 2, tol = 1e-10; end
[V,D] = eig(A);
lam   = diag(D);

% Identify the real eigenvalue
isReal = abs(imag(lam)) <= tol;
if sum(isReal) ~= 1
    error('Expected exactly one real eigenvalue for a 3x3 with one real + one complex pair.');
end
iReal = find(isReal,1);

% Pick the complex eigenvalue with positive imag part
idxCmp = setdiff(1:3, iReal);
if imag(lam(idxCmp(1))) > 0
    iCmp = idxCmp(1);
else
    iCmp = idxCmp(2);
end

q1  = real(V(:,iReal));
vc  = V(:,iCmp);
u2  = real(vc);
v2  = imag(vc);

% Normalize and Gram-Schmidt clean-up for the plane axes
q1  = nz_norm(q1, tol);
u2  = nz_norm(u2, tol);
v2  = v2 - (u2' * v2) * u2;
v2  = nz_norm(v2, tol);

Q = [q1, u2, v2];
if abs(det(Q)) < 1e-12
    error('Modal basis Q is near-singular.');
end
U2 = Q(:,2:3);
end

function x = nz_norm(x, tol)
if nargin<2, tol = 1e-10; end
n = norm(x);
if n >= tol
    x = x / n;
end
end

function fig = make_figure_for_system(A, Q, q1, U2, x0_a, x0_b, tag, ...
                                      T, N, savepath, doShow, span, planeAlpha, useEig)
t = linspace(0, T, N).';

% Projections (row 2)
pL1_a = proj_on_line(x0_a, q1);
pL2_a = proj_on_plane(x0_a, U2);
pL1_b = proj_on_line(x0_b, q1);
pL2_b = proj_on_plane(x0_b, U2);

% Trajectories (row 3)
Xa = trajectory(A, x0_a, t, useEig);
Xb = trajectory(A, x0_b, t, useEig);

visFlag = 'on'; if ~doShow, visFlag = 'off'; end
fig = figure('Color','w','Visible',visFlag,'Renderer','opengl');
tl  = tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');

% Row 1: subspaces only
ax1 = nexttile(tl, [1 2]); hold(ax1,'on'); grid(ax1,'on');
draw_subspaces(ax1, q1, U2, planeAlpha, span);
title(ax1, sprintf('%s: invariant subspaces (L_1 line, L_2 plane)', tag), 'Interpreter','tex');
xlabel(ax1,'x_1'); ylabel(ax1,'x_2'); zlabel(ax1,'x_3');
equalize3d(ax1);

% Row 2 left: projections for x0^(a)
ax21 = nexttile(tl, 3); hold(ax21,'on'); grid(ax21,'on');
draw_subspaces(ax21, q1, U2, 0.18, span);
scatter3(ax21, x0_a(1), x0_a(2), x0_a(3), 35, 'filled', 'DisplayName','x_0^{(a)}');
scatter3(ax21, pL1_a(1), pL1_a(2), pL1_a(3), 25, 'filled', 'DisplayName','proj L_1');
scatter3(ax21, pL2_a(1), pL2_a(2), pL2_a(3), 25, 'filled', 'DisplayName','proj L_2');
plot3(ax21, [x0_a(1) pL1_a(1)], [x0_a(2) pL1_a(2)], [x0_a(3) pL1_a(3)], '--', 'LineWidth',1.2);
plot3(ax21, [x0_a(1) pL2_a(1)], [x0_a(2) pL2_a(2)], [x0_a(3) pL2_a(3)], '--', 'LineWidth',1.2);
title(ax21, sprintf('%s: projections of x_0^{(a)} onto L_1 and L_2', tag), 'Interpreter','tex');
legend(ax21,'Location','northwest','Box','off');
equalize3d(ax21);

% Row 2 right: projections for x0^(b)
ax22 = nexttile(tl, 4); hold(ax22,'on'); grid(ax22,'on');
draw_subspaces(ax22, q1, U2, 0.18, span);
scatter3(ax22, x0_b(1), x0_b(2), x0_b(3), 35, 'filled', 'DisplayName','x_0^{(b)}');
scatter3(ax22, pL1_b(1), pL1_b(2), pL1_b(3), 25, 'filled', 'DisplayName','proj L_1');
scatter3(ax22, pL2_b(1), pL2_b(2), pL2_b(3), 25, 'filled', 'DisplayName','proj L_2');
plot3(ax22, [x0_b(1) pL1_b(1)], [x0_b(2) pL1_b(2)], [x0_b(3) pL1_b(3)], '--', 'LineWidth',1.2);
plot3(ax22, [x0_b(1) pL2_b(1)], [x0_b(2) pL2_b(2)], [x0_b(3) pL2_b(3)], '--', 'LineWidth',1.2);
title(ax22, sprintf('%s: projections of x_0^{(b)} onto L_1 and L_2', tag), 'Interpreter','tex');
legend(ax22,'Location','northwest','Box','off');
equalize3d(ax22);

% Row 3 left: trajectory from x0^(a)
ax31 = nexttile(tl, 5); hold(ax31,'on'); grid(ax31,'on');
draw_subspaces(ax31, q1, U2, 0.12, span);
plot3(ax31, Xa(:,1), Xa(:,2), Xa(:,3), 'LineWidth',2.0, 'DisplayName','x(t) from x_0^{(a)}');
scatter3(ax31, x0_a(1), x0_a(2), x0_a(3), 30, 'filled');
title(ax31, sprintf('%s: trajectory from x_0^{(a)}', tag), 'Interpreter','tex');
xlabel(ax31,'x_1'); ylabel(ax31,'x_2'); zlabel(ax31,'x_3');
legend(ax31,'Box','off');
equalize3d(ax31);

% Row 3 right: trajectory from x0^(b)
ax32 = nexttile(tl, 6); hold(ax32,'on'); grid(ax32,'on');
draw_subspaces(ax32, q1, U2, 0.12, span);
plot3(ax32, Xb(:,1), Xb(:,2), Xb(:,3), 'LineWidth',2.0, 'DisplayName','x(t) from x_0^{(b)}');
scatter3(ax32, x0_b(1), x0_b(2), x0_b(3), 30, 'filled');
title(ax32, sprintf('%s: trajectory from x_0^{(b)}', tag), 'Interpreter','tex');
xlabel(ax32,'x_1'); ylabel(ax32,'x_2'); zlabel(ax32,'x_3');
legend(ax32,'Box','off');
equalize3d(ax32);

sgtitle(tl, sprintf('%s — Subspaces, projections, and trajectories', tag), 'FontWeight','bold');

% Save
if ~isempty(savepath)
    try
        exportgraphics(fig, savepath, 'Resolution',180);
    catch
        % Fallback if exportgraphics unavailable
        saveas(fig, savepath);
    end
end

end

function X = trajectory(A, x0, t, useEig)
% Return X(t) as [N x 3], evaluating e^{At} x0
N = numel(t);
X = zeros(N,3);
if nargin<4, useEig = false; end
if useEig
    [V,D] = eig(A);
    Vinv  = inv(V);
    for i = 1:N
        M = V * diag(exp(diag(D)*t(i))) * Vinv;
        X(i,:) = real(M * x0).';
    end
else
    for i = 1:N
        X(i,:) = (expm(A * t(i)) * x0).';
    end
end
end

function p = proj_on_line(x, a)
% Orthogonal projection of x onto span(a)
c = (x.'*a) / (a.'*a);
p = c * a;
end

function p = proj_on_plane(x, U)
% Orthogonal projection of x onto span(U), U in R^{3x2}
coef = (U.'*U) \ (U.'*x);
p = U * coef;
end

function draw_subspaces(ax, q1, U2, alphaPlane, span)
% L2 plane
a = linspace(-span, span, 10);
b = linspace(-span, span, 10);
[Aa,Bb] = meshgrid(a,b);
P = U2(:,1)*Aa + U2(:,2)*Bb;  % 3 x m x n via implicit expansion
surf(ax, P(1,:), P(2,:), P(3,:), ...
     'EdgeColor','none', 'FaceAlpha',alphaPlane);
text3(ax, 0.9*(U2(1,1)+U2(1,2)), 0.9*(U2(2,1)+U2(2,2)), 0.9*(U2(3,1)+U2(3,2)), ...
      'L_2 plane', 'Interpreter','tex');

% L1 line
pts = [-span, span];
L1  = q1.' .* pts.';
plot3(ax, L1(:,1), L1(:,2), L1(:,3), '--', 'LineWidth',2.0);
text3(ax, 1.05*q1(1), 1.05*q1(2), 1.05*q1(3), 'L_1 line', 'Interpreter','tex');

axis(ax,'vis3d'); box(ax,'on');
end

function equalize3d(ax)
% Make a cubic box that encloses current xyz limits
xl = xlim(ax); yl = ylim(ax); zl = zlim(ax);
mins = [xl(1) yl(1) zl(1)];
maxs = [xl(2) yl(2) zl(2)];
ctr = (mins + maxs)/2;
rad = max(maxs - mins)/2;
xlim(ax, ctr(1)+rad*[-1 1]);
ylim(ax, ctr(2)+rad*[-1 1]);
zlim(ax, ctr(3)+rad*[-1 1]);
end
