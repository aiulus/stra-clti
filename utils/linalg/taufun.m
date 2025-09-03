function tau = taufun(Y, S)
%TAUFUN Smoothed condition number variant: cond(Y*S*Y').
tau = cond(Y*S*Y.');
end
