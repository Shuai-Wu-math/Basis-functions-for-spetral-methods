function [dy,y] = lpsm(n,x)
% This code computes Legendre polynomials and their derivatives.
% It takes two inputs: n, which specifies the order up to which the polynomials 
% are calculated, and x, the specific point at which to evaluate them. 
% If the output quantity is 1, it returns the polynomial value denoted as y. 
% If the output quantity is 2, it returns both the polynomial value y and 
% the derivative value dy.

if nargout==1
    
     if n==0, y=ones(size(x));  return; end
     if n==1, y=[ones(size(x));x]; return; end
     
     poly1=ones(size(x)); poly2=x;   
     y=[poly1;poly2];
     for  k=2:n                      % Three-term recurrence relation:  
	   polyn=((2*k-1)*x.*poly2-(k-1)*poly1)/k; 
       poly1=poly2; poly2=polyn; y=[y;polyn];	
     end
     
end

if nargout==2
     if n==0, y=ones(size(x)); dy=zeros(size(x)); return;end
     if n==1, y=[ones(size(x));x]; 
              dy=[zeros(size(x));ones(size(x))]; return; end

     poly1=ones(size(x)); pderl=zeros(size(x));poly2=x; pder2=ones(size(x));
     y=[poly1;poly2]; dy=[pderl;pder2];
    for k=2:n                          
      polyn=((2*k-1)*x.*poly2-(k-1)*poly1)/k; 
      pdern=pderl+(2*k-1)*poly2;  
 	  poly1=poly2; poly2=polyn; y=[y;polyn];
	  pderl=pder2; pder2=pdern; dy=[dy;pdern];
    end
      
end

end