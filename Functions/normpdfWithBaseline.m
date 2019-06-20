function y = normpdfWithBaseline(x,params)
%


A = params(1);
mu = params(2);
sigma = params(3);
C=params(4);
y= A*exp(-0.5.*((x-mu)./(sigma)).^2) +C;
