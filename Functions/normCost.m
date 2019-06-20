function fVal = normCost(x,A,mu,sig,C)

params =[A mu sig C];
y=normpdfWithBaseline(x,params);

fVal = -sum(log(y));
