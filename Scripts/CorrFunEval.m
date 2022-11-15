function f = CorrFunEval(type, x, y, a)
%		   
%		  f = CorrFunEval(type, x, y, a)
%		  
%		  Compute the autocorrelation function 
%			
%		  type  : 'exp' for exponential, 'exp2' for exponential-square
%		  x , y : points of evaluation
%		  a     : length of correlation
%		  
%		  NB    : this function is vectorized, i.e x,y,a can be vectors 
  
  switch type
   case 'exp'
	f = exp(- sum(abs(x-y)./a)) ;
   case 'exp2'
	f = exp(- sum(((x-y).^2)./ (a.^2))) ;
   otherwise 
	fprintf('\n This correlation function type is not supported \n');
	error(' ' );
  end;
  
  
