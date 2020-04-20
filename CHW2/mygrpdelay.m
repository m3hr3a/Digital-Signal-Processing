function [w,nd] = mygrpdelay(h)
N = length(h) * 100 ; 
H = (fft(h,N+1)) ; 
w = 0 : 2*pi/(N+1) : 2*pi - 2*pi/(N+1) ; 
dw = w(2:end) - w(1 : end-1);
dH = H(2:end) - H(1 : end-1);
nd = 0 - real(1i *dH./dw ./ H(1:end-1));
w = w(1:end-1);
