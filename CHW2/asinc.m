function [y] = asinc( L , w )
 y = sin(1/2*w*L)./sin(w/2);
 kmax = ceil(max(w)/2/pi); 
 kmin = floor(max(w)/2/pi); 
 for k = kmin : 1 : kmax % this loop check zeros div.
   zerodiv = find(abs(w-2*pi*k)<1e-4);
   % correct zero div. s with lhospital theorem 
   y(zerodiv) = L * cos(w(zerodiv)*L/2) ./ cos(w(zerodiv)/2);
 end
   