function [xi,wi] = lpsnw(n)
% This function computes the Legendre-Gauss quadrature nodes and weights.
% It takes an input n and returns the zeros of the nth-degree Legendre polynomial 
% along with their corresponding weights. The n quadrature nodes ensure an 
% algebraic precision of 2n+1.

thetak=(4*(1:n)-1)*pi/(4*n+2);
xi=-(1-(n-1)/(8*n^3)-(39-28./sin(thetak).^2)/(384*n^4)).*cos(thetak);
ep=eps*10;                            % error tolerance for stopping iteration
ze1=xi+ep+1;
 
while max(abs(ze1-xi))>=ep            % Newton's iteration procedure
    ze1=xi;
    [dy,y]=lps(n,xi);
    xi=xi-y./dy;  
end                                  

wi=(2./((1-xi.^2).*dy.^2))';

end

