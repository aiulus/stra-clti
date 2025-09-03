function makeFigure_Ax0(A_sys, q1, U2, x0_a, x0_b, tag, savepath, T, N)
%MAKEFIGURE_AX0 Build the 3x2 layout and draw everything.

if nargin < 9 || isempty(T), T = 6.0; end
if nargin < 10 || isempty(N), N = 600; end

t = linspace(0, T, N);

% --- Projections (row 2)
pL1_a = projOnLine(x0_a, q1);
pL2_a = projOnPlane(x0_a, U2);

pL1_b = projOnLine(x0_b, q1);
pL2_b = projOnPlane(x0_b, U2);

% --- Trajectories (row 3)
Xa = grid_trajectory(A_sys, x0_a, t);
Xb = grid_trajectory(A_sys, x0_b, t);

% --- Layout: tiledlayout 3x2; row 1 spans two columns
fig = figure('Color','w','Position',[100 100 1100 1000]);
tl = tiledlayout(fig, 3, 2, 'TileSpacing','compact', 'Padding','compact');

% Row 1: subspaces only
ax1 = nexttile(tl, [1 2]); hold(ax1,'on'); view(ax1,3);
drawSubspaces(ax1, q1, U2, 0.22, 2.4);
title(ax1, sprintf('%s: invariant subspaces (L_1 line, L_2 plane)', tag));
xlabel(ax1,'x_1'); ylabel(ax1,'x_2'); zlabel(ax1,'x_3');
equal3d(ax1); grid(ax1,'on');

% Row 2, col 1: projections for x0^(a)
ax21 = nexttile(tl, 3); hold(ax21,'on'); view(ax21,3);
drawSubspaces(ax21, q1, U2, 0.18, 2.4);
scatter3(ax21, x0_a(1), x0_a(2), x0_a(3), 35, 'filled', 'DisplayName','x_0^{(a)}');
scatter3(ax21, pL1_a(1), pL1_a(2), pL1_a(3), 25, 'filled', 'MarkerFaceAlpha',0.8);
scatter3(ax21, pL2_a(1), pL2_a(2), pL2_a(3), 25, 'filled', 'MarkerFaceAlpha',0.8);
plot3(ax21, [x0_a(1) pL1_a(1)], [x0_a(2) pL1_a(2)], [x0_a(3) pL1_a(3)], 'k--', 'LineWidth',1.1);
plot3(ax21, [x0_a(1) pL2_a(1)], [x0_a(2) pL2_a(2)], [x0_a(3) pL2_a(3)], 'k--', 'LineWidth',1.1);
title(ax21, sprintf('%s: projections of x_0^{(a)} onto L_1 and L_2', tag));
legend(ax21,'Location','northwest','Box','off'); grid(ax21,'on'); equal3d(ax21);
xlabel(ax21,'x_1'); ylabel(ax21,'x_2'); zlabel(ax21,'x_3');

% Row 2, col 2: projections for x0^(b)
ax22 = nexttile(tl, 4); hold(ax22,'on'); view(ax22,3);
drawSubspaces(ax22, q1, U2, 0.18, 2.4);
scatter3(ax22, x0_b(1), x0_b(2), x0_b(3), 35, 'filled', 'DisplayName','x_0^{(b)}');
scatter3(ax22, pL1_b(1), pL1_b(2), pL1_b(3), 25, 'filled', 'MarkerFaceAlpha',0.8);
scatter3(ax22, pL2_b(1), pL2_b(2), pL2_b(3), 25, 'filled', 'MarkerFaceAlpha',0.8);
plot3(ax22, [x0_b(1) pL1_b(1)], [x0_b(2) pL1_b(2)], [x0_b(3) pL1_b(3)], 'k--', 'LineWidth',1.1);
plot3(ax22, [x0_b(1) pL2_b(1)], [x0_b(2) pL2_b(2)], [x0_b(3) pL2_b(3)], 'k--', 'LineWidth',1.1);
title(ax22, sprintf('%s: projections of x_0^{(b)} onto L_1 and L_2', tag));
legend(ax22,'Location','northwest','Box','off'); grid(ax22,'on'); equal3d(ax22);
xlabel(ax22,'x_1'); ylabel(ax22,'x_2'); zlabel(ax22,'x_3');

% Row 3, col 1: trajectory from x0^(a)
ax31 = nexttile(tl, 5); hold(ax31,'on'); view(ax31,3);
drawSubspaces(ax31, q1, U2, 0.12, 2.4);
plot3(ax31, Xa(:,1), Xa(:,2), Xa(:,3), 'LineWidth',2.0, 'DisplayName','x(t) from x_0^{(a)}');
scatter3(ax31, x0_a(1), x0_a(2), x0_a(3), 30, 'filled');
title(ax31, sprintf('%s: trajectory from x_0^{(a)}', tag));
legend(ax31,'Location','northwest','Box','off'); grid(ax31,'on'); equal3d(ax31);
xlabel(ax31,'x_1'); ylabel(ax31,'x_2'); zlabel(ax31,'x_3');

% Row 3, col 2: trajectory from x0^(b)
ax32 = nexttile(tl, 6); hold(ax32,'on'); view(ax32,3);
drawSubspaces(ax32, q1, U2, 0.12, 2.4);
plot3(ax32, Xb(:,1), Xb(:,2), Xb(:,3), 'LineWidth',2.0, 'DisplayName','x(t) from x_0^{(b)}');
scatter3(ax32, x0_b(1), x0_b(2), x0_b(3), 30, 'filled');
title(ax32, sprintf('%s: trajectory from x_0^{(b)}', tag));
legend(ax32,'Location','northwest','Box','off'); grid(ax32,'on'); equal3d(ax32);
xlabel(ax32,'x_1'); ylabel(ax32,'x_2'); zlabel(ax32,'x_3');

sgtitle(tl, sprintf('%s — Subspaces, projections, and trajectories', tag), 'FontWeight','bold');

% Save & show
if ~isempty(savepath)
    try
        exportgraphics(fig, savepath, 'Resolution', 180);
    catch
        saveas(fig, savepath);
    end
end

end
