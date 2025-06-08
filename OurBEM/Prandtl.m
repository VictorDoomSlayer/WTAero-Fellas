function [aPr,aprimePr,ftotal] = Prandtl(B,r,R,lambda,a,aprime,rRoot)

mu = r/R;
muRoot = rRoot/R;
ftip = (2/pi)*acos(exp(-0.5*B*((1-mu)/mu)*sqrt(1+(((lambda^2) * mu^2)/(1-a)^2))));
froot = (2/pi)*acos(exp(-0.5*B*((mu-muRoot)/mu)*sqrt(1+(((lambda^2) * mu^2)/(1-a)^2))));
ftotal = ftip*froot;

aPr = a/ftotal;
aprimePr = aprime/ftotal;
end