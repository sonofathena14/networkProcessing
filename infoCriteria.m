% Input: k = # of parameters
%        n = number of data points
%        J = cost (SSE)
% Output: AIC and BIC scores

function [AIC, BIC] = infoCriteria(n, k, J)

AIC = 2*k + 2*log(J);
BIC = log(n)*k + 2*log(J);

end