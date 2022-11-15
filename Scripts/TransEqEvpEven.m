function z = TransEqEvpEven(x,c,a)
%
%     y = TransEqEvpEven(x,c,a)
%     
%     To solve the eigenvalue problem associated with the exponential
%     kernel, some solutions of transcendant algebraic equations are needed
%     
%     See book by Spanos and Ghanem, 1991
%     
   z = c * tan (a*x) + x;
return
