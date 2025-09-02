function J = Jordan(eigvals)
% Constructs real Jordan blocks from a list of eigenvalues
    eigv = eigvals(:);
    eigv = unique([eigv; conj(eigv)], 'stable'); % ensure conjugate closure

    eigv = eigv(imag(eigv) >= 0); % remove eigenvalues with negative imaginary parts
    ndim = 0;
    for k = 1:numel(eigv)
        if abs(imag(eigv(k))) < eps
            ndim = ndim + 1;
        else
            ndim = ndim + 2;
        end
    end
    J = zeros(ndim);
    dj = 1;
    for k = numel(eigv)
        lambda_k = eigv(k);
        if abs(imag(lambda_k)) < eps
            J(dj, dj) = real(lambda_k);
            dj = dj + 1;
        else
            a = real(lambda_k);
            b = imag(lambda_k);
            J(dj:dj+1, dj:dj+1) = [a, b; -b, a];
            dj = dj + 2;
        end
    end
    J = real(J);
end