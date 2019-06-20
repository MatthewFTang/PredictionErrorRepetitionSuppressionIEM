function y = normpdfWithBaseline1(x,A,mu,sigma,C)
%


y= A*exp(-0.5.*((x-mu)./(sigma)).^2) +C;
