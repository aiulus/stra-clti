function p = projOnLine(x, a)
%PROJONLINE Orthogonal projection of x onto span(a).
c = (a(:)'*x(:)) / (a(:)'*a(:));
p = c * a(:);
end
