function RF  =  KLExpansion(SizeDomain,CorrLength,OrderExp)
%
%      RF  =  KLExpansion(SizeDomain,CorrLength,OrderExp)
%
%  Compute the eigenvalues and eigenfunction parameters 
%  for the problem associated with the Karhunen-Loeve 
%  expansion of the exponential kernel.
%  (1D and 2D assuming same correlation length in both directions)
%
%  ----Input----
%
%          a  : Domain of definition of the 1D eigenvalue problem is [-a : a]
%                 !!!! This may not be the domain on which the random field
%		          is discretized. To limit variance error at the boundaries 
%			      indeed, a larger domain may be defined.
%
%  CorrLength :   correlation length of the process, that is
%                 the parameter in the correlation function
%                 (c = 1/ CorrLength)
%
%
%    NbEig1D  : number of eigenvalues calculated for 1D problem
%    NbEig2D  : number of eigenvalues calculated for 2D problem
%  
%  ----Output----
%
%  CorrEVP    : is a [NbEig1D , 3] matrix. Each row n contains the three
%               following coefficients: 
%
%	            wn : solution of the transcendental equation in w
%			        c - w * tan(w*a) = 0   if n odd
%			        c * tan(w*a) +w  = 0   if n even
%
%	            lambda : 2*c/( wn^2 +c^2) is the n-th eigenvalue 
%  
%               alpha: the normalization factor for the eigenfunction
%       			 (the latter is cos (wn*x)  if n odd
%                     sin (wn*x)  if n even)   
%
%  Eigs2DValue: contains the NbEig2D eigenvalues obtained by multiplying
%               any two 1D eigenvalues. Ordered in descending order
%
%  Eigs2Dpos  : for each 2D eigenvalue, gives both 1D indices.
%




%--------------------------------------------------------------------  
%       Initializations
%--------------------------------------------------------------------
optimset.Display = 'off';
opti = optimset;
CorrEVP = zeros([OrderExp , 3]);  

dimension = length(SizeDomain);   % Problem dimension

switch dimension
   
case 1     % KL expansion in one dimension
   a = SizeDomain(1);
   c = 1 / CorrLength (1);
   %--------------------------------------------------------------------  
   %       Generating the matrix CorrEVP
   %--------------------------------------------------------------------
   
   for i = 0 : ceil(OrderExp/2) 
      % compute interval where zeros are to be found
      intv = [max((2*i-1)*pi/(2*a)+0.00000001, 0) (2*i+1)*pi/(2*a)-0.00000001 ];
      
      % compute data associated with equation : c * tan (a*x) + x
      if ((i > 0) & (2*i <=OrderExp))
         
            % Matlab 5.3 version
            wnstar = fzero('TransEqEvpEven' , intv,opti,c,a);
            % Matlab 5.0 - 5.2 version
            %wnstar = fzero('TransEqEvpEven' , intv,[],[],c,a);
            
            CorrEVP(2 * i , 1) = wnstar;                                  % omega_n
            CorrEVP(2 * i , 2) = 2*c/( wnstar^2 +c^2);                    % lambda_n
            CorrEVP(2 * i , 3) = 1/sqrt(a - sin(2*wnstar*a)/(2*wnstar));  % alpha_n
         end;
         
         % compute data associated with equation : c - x * tan (a*x) = 0
         if ((2*i +1) <= OrderExp)
            
            % Matlab 5.3 version
            wn = fzero('TransEqEvpOdd' , intv ,opti,c,a);
            % Matlab 5.0 - 5.2 version
            %wn = fzero('TransEqEvpOdd' , intv ,[],[],c,a);
            
            CorrEVP(2 * i +1 , 1) = wn ;                              % omega_n
            CorrEVP(2 * i +1 , 2) = 2*c/( wn^2 +c^2);                 % lambda_n
            CorrEVP(2 * i +1 , 3) = 1/sqrt(a + sin(2*wn*a)/(2*wn));   % alpha_n
         end;
      end;
      RF.EVPx = CorrEVP;
      RF.Eigs = CorrEVP(: , 2);
      
   case 2
      %------------------------------------------------------------------------
      %        Computing the 1D basis data for different domain size and
      %        correlation 
      %------------------------------------------------------------------------
      
      ax = SizeDomain(1); ay = SizeDomain(2);
      corrx = CorrLength (1); corry = CorrLength (2);
      RF.EVPx = getfield(KLExpansion(ax,corrx,OrderExp) , ...
         'EVPx');
      if ((abs(ax -ay) < 1e-10) & ( abs (corrx - corry) < 1e-10))
         RF.EVPy = RF.EVPx ;
      else
         RF.EVPy = getfield(KLExpansion(ay,corry,OrderExp) ,'EVPx');
      end;
      
      
      %------------------------------------------------------------------------
      %        Computing the 2D basis
      %------------------------------------------------------------------------
      %     2D problem is defined over [-ax ; ax ] * [-ay ; ay]
      %     Correlation length are 
      %
      %     2D eigenvalues are product of 1D eigenvalues. They are ranked in
      %     descending order.
      %     Each 2D egv is designed by the index [i,j] of the product
      %     This is stored in a cell array of strutures 
      
      
      Eigs1Dx = RF.EVPx(:, 2);                   % 1D eigenvalues
      Eigs1Dy = RF.EVPy(:, 2);                   % 1D eigenvalues
      
      EigsProduct = (Eigs1Dx * Eigs1Dy');        % 2D eigenvalues with double
      % counting 
      
      
      for i = 1 : OrderExp
         [MaxColumn , IndRow] = max(EigsProduct);   
         % MaxColumn is an array, each term
         % being the max of EigsProduct(:,j)
         % IndRow(j) gives the number of
         % the row wher the max is
         % encountered
         [eigmax	, jj ] = max(MaxColumn);    % eigmax is the largest eigenvalue
         % remaining in Eigs2D 
         %    if (jj > MaxIndColumn) 
         %	  MaxIndColumn = jj;
         %	end;
         ii = IndRow(jj);         									
         Eigs2Dvalue(i) = eigmax;
         Eigs2Dpos(i, :) = [ii,jj];
         EigsProduct(ii,jj) = 0;          % set to zero the current max to find
         % the following
      end;
      
      
      RF.Eigs = Eigs2Dvalue;
      RF.EVP2pos = Eigs2Dpos;
      
      
   otherwise
   end;
   
   return
   
 