function[dc_fac_H,xn] = filter_abaqus_heaviside(beta,x)
xn = zeros(size(x));
dc_fac_H = zeros(size(x));
%diff = 1;
l2 = 1;
l1 = 0;
while (l2-l1) > 0.0001
eta = (l2+l1)/2;
k = find(x<=eta);
k1 = find(x>eta);
xn(k) = eta*(exp(-beta*(1-x(k)/eta))-(1-x(k)/eta)*exp(-beta));
xn(k1) = (1-eta)*(1-exp(-beta*(x(k1)-eta)/(1-eta))+(x(k1)-eta)*exp(-beta)/(1-eta))+eta;
diff = sum(xn)-sum(x);
if diff == 0
    l2 =eta;
    l1 = eta;
elseif diff > 0
    l1 = eta;
else
    l2 = eta;
end
end
dc_fac_H(k) = beta*exp(-beta*(1-x(k)/eta))+exp(-beta);
dc_fac_H(k1) = beta*exp(-beta*((x(k1)-eta)/(1-eta)))+exp(-beta);
% may 22, 2020
dc_fac_H = repmat(dc_fac_H',length(x),1);