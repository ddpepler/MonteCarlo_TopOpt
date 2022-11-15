function y = TransEqEvpOdd(x,c,a)
%
%     y = TransEqEvpOdd(x,c,a)
%     
%     To solve the eigenvalue problem associated with the exponential
%     kernel, some solutions of transcendant algebraic equations are needed
%     
%     See book by Spanos and Ghanem, 1991
%     
  y = c - x * tan (a*x);
return
