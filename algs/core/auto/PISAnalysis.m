function out = PISAnalysis(Y, S, L)
%PISANALYSIS Bundle practical identifiability scores (SCN, PIS, kappa).

out = struct();
out.PIS   = wstarfun(Y,S,L);
out.kappa = kappafun(Y);
out.SCN   = taufun(Y,S);
end
