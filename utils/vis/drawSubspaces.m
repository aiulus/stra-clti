function drawSubspaces(ax, q1, U2, alphaPlane, span)
%DRAWSUBSPACES Render L2 plane (span of U2) and L1 line (span of q1).
if nargin < 4, alphaPlane = 0.22; end
if nargin < 5, span = 2.4; end

% L2 plane as parametric surface: a*U2(:,1) + b*U2(:,2)
aVals = linspace(-span, span, 16);
bVals = linspace(-span, span, 16);
[AA, BB] = meshgrid(aVals, bVals);
Xpl = U2(1,1)*AA + U2(1,2)*BB;
Ypl = U2(2,1)*AA + U2(2,2)*BB;
Zpl = U2(3,1)*AA + U2(3,2)*BB;
s = surf(ax, Xpl, Ypl, Zpl, 'EdgeColor','none', 'FaceAlpha', alphaPlane);
uistack(s,'bottom');

% label L2
txtL2 = 0.9*U2(:,1) + 0.9*U2(:,2);
text(ax, txtL2(1), txtL2(2), txtL2(3), 'L_2 plane');

% L1 line
pts = [-span; span];
L1  = [pts*q1(1), pts*q1(2), pts*q1(3)];
plot3(ax, L1(:,1), L1(:,2), L1(:,3), 'k--', 'LineWidth', 2.0);
text(ax, 1.05*q1(1), 1.05*q1(2), 1.05*q1(3), 'L_1 line');
end