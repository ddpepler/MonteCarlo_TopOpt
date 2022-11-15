function f = EvalBasisFunction(RF,X,n)
%
%  f =   EvalBasisFunction(RF,X,n)
% 
%  A discretized random filed is represented by :
%  H_app(x) = mu + Sum_{i=1..OrderExp)  h_j(x) xi_j
%  
%	 where mu       : mean value of homogeneous field
%		   OrderExp : number of terms in series expansion
%		   xi_j     : uncorrelated N(0,1) variates
%		   h_j(x)   : some deterministic basis function depending on ...
%					   the discretization scheme.
%
%  this function evaluate the n-th basis function h_n at point X:
%		   RF       : the random field structure, obtained by DiscRandomField
%		   X        : the current point
%		   n        : the basis function ID
%  
%			 

  switch RF.DiscScheme
   case 'KL'
	% The current point X is defined in user coordinates.
	% EvalPhi works in a basis so that the structure
	% is enclosed in a rectangle centered on the origin
	% The relationship between theses basis is a translation T. 
	%
	%        sqrt(RF.Eigs(n)) is the sqrt of corresponding eigenvalue
	%        RF.Stdv is the Std deviation
	
	f = RF.Stdv * sqrt(RF.Eigs(n)) * ...
		EvalPhi(RF,X-RF.Translation,n);
	
	
   case 'EOLE'
	LL = length(RF.COORD);
	Rho_vV =   zeros([LL , 1]);   % Vector containing correlation between
                                 % current point X and Xi= nodes of the
                                 % RF mesh.
	for i = 1 : LL
	  Rho_vV(i) = CorrFunEval(RF.CorrType, X , RF.COORD(i,:) , ...
							  RF.CorrLength);
	end;
	Phi_n = RF.Phi(:,n);
	%Rho_vV	
	f = RF.Stdv * (Phi_n' * Rho_vV) / sqrt(RF.Eigs(n));
	
   case 'OSE'
	%   Compute h_i(x) = sqrt((2*i+1)/(2*a)) P_i((x-Y)/a)   i =0, ...
	%   !!!! P_0 is actually stored in Hn(1)
	%
	Hn = OSEBasis(RF,X);
	f = 0;
	
	for i =1 : RF.OrderExp
	  f = f +  RF.Phi(i,n) * Hn(i);
	end
	f =f * RF.Stdv * sqrt(RF.Eigs(n));
   otherwise
  end;
  
  return
  
  
  
  
  
%----------------------------------------------------------------------------  
  function f = EvalPhi(RFData,X,n)
%
%     	f = EvalPhi(RFData,X,n) 
%     	
%     	Evaluate the n-th eigenfunction of the correlation kernel at point X
%     	in 1D or 2D
%
%       !!!!!!!! f does not encompass the stdv factor, neither the eigenvalue
%
%       RFData : Data structure containing the deterministic basis
%       description. So far allows only Karhunen Loeve expansion.
%       in 1D  : RFData has the sole field 'EVPx'
%       in 2D  : RFData has 4 fields 'EVPx', 'EVPy', Eigs, EVP2pos
%
%     	X    : either scalar (-->x ) or vector of length 2 (-->[x,y])
%
%     	n    : order of the eigenfunction
  
  
switch length(X)
 case 1
  EVP1     = RFData.EVPx   ;
  w      = EVP1(n,1);
  lambda = EVP1(n,2);
  alpha  = EVP1(n,3);
  if (mod(n,2) == 1)    % odd rank eigenfunction
	f = alpha * cos(w*X);
  else                 % even rank eigenfunction
	f = alpha * sin(w*X);
  end;

  
 case 2
  EVPx     = RFData.EVPx   ;
  EVPy     = RFData.EVPy   ;
  Eigs  = RFData.Eigs;
  EVP2pos  = RFData.EVP2pos; 

  x = X(1);            % extracting coordinates
  y = X(2);
  ii = EVP2pos(n,1);
  jj = EVP2pos(n,2);
  
  %  Direct product by recursive call
  % f =  EvalBasisFunction(x,ii) * EvalBasisFunction(y,jj) ;
  
  %  Direct product
  w      = EVPx(ii,1);
  lambda = EVPx(ii,2);
  alpha  = EVPx(ii,3);
  if (mod(ii,2) == 1)    % odd rank eigenfunction
  		f = alpha * cos(w*x);
  else                 % even rank eigenfunction
  		f = alpha * sin(w*x);
  end;
  
  w      = EVPy(jj,1);
  lambda = EVPy(jj,2);
  alpha  = EVPy(jj,3);
  if (mod(jj,2) == 1)    % odd rank eigenfunction
  		f =  f * alpha * cos(w*y);
  else                 % even rank eigenfunction
  		f =  f * alpha * sin(w*y);
  end;

 otherwise 
  fprintf('Error in calling EvalBasisFunction \n');
  fprintf('  current point x has 1 or 2 coordinates');
end

return
%---------------------------------------------------------------------------
