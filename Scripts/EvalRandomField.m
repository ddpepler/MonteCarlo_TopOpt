function f = EvalRandomField(RF,X,Xi)
%  
%  f = EvalRandomField(RF,X,Xi)
%  
%  Evaluate a realization of a discretized random field at point X.
%  
%  RF    : the random field structure, obtained by DiscRandomField
%  X     : the current point
%  Xi    : the vector of outcome of N(0,1) appearing in the expansion : 
%			H_app(x) = mu + Sum_{i=1..OrderExp)  h_j(x) Xi_j
%
%  NB : if the field is lognormal, this function returns
%         exp(mu + Sum_{i=1..OrderExp)  h_j(x) Xi_j)


f = RF.Mean;
for i = 1 : length(Xi)
   f = f + Xi(i) *  EvalBasisFunction(RF,X,i);
end;

switch RF.Type
case 'Gaussian'
   % Nothing
case 'Lognormal'
   f = exp(f) ;
end
