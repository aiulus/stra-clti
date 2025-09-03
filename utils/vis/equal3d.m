function equal3d(ax)
%EQUAL3D Make 3D axes use equal data scaling.
% Works across MATLAB versions better than axis equal alone.
xl = xlim(ax); yl = ylim(ax); zl = zlim(ax);
c = [mean(xl) mean(yl) mean(zl)];
r = max([diff(xl) diff(yl) diff(zl)]) / 2;
xlim(ax, c(1)+[-r r]); ylim(ax, c(2)+[-r r]); zlim(ax, c(3)+[-r r]);
pbaspect(ax,[1 1 1]);
end
