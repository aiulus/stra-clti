function p = projOnPlane(x, U)
%PROJONPLANE Orthogonal projection of x onto span(U) (U is n x 2).
G = U.'*U;
coef = G \ (U.'*x(:));
p = U*coef;
end