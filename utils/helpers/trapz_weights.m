function w = trapz_weights(t)
n = numel(t);
w = zeros(1,n);
w(1)   = 0.5*(t(2)-t(1));
w(end) = 0.5*(t(end)-t(end-1));
for k = 2:n-1
    w(k) = 0.5*(t(k+1)-t(k-1));
end
end