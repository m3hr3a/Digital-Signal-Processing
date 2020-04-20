%% DSP - HW3 - Programming part -1st MATLAB assignment
%% Instructor : Dr. Babaie-zadeh 
%% Student : Mehrsa Pourya 
%% project 1
%% Ex. 1.1
clear 
clc
close all
% a :
a = [ 1 0 0.9 ]  % y coefs
b = [ 0.3 0.6 0.3 ]  % x coefs
% b : 
analres=zeros(1, 127); % allocate analytic response
% define analytic response based on calculations report in report file
for n = 1 : 127 
    if( mod(n, 4)==0)
        analres(n) = -1/30 * (sqrt(0.9))^n ; 
    end
    if( mod(n, 4)==1)
        analres(n) = 0.6 * (sqrt(0.9))^(n-1) ; 
    end
    if( mod(n, 4)==2)
        analres(n) = 1/30 * (sqrt(0.9))^n ; 
    end
    if( mod(n, 4)==3)
        analres(n) = -0.6 * (sqrt(0.9))^(n-1) ; 
    end
end
analres = [0.3 analres]; 
% analytic impulse response stem
subplot(2,1,1) 
stem(0:127,analres)
grid on 
xlabel('n') 
ylabel('y[n]')
title('Analytical impulse response visualization')
% c : simulated impulse response
imp = [ 1 zeros(1, 127) ] ;
impresponse = filter(b ,a , imp) ;   % calculate impulse res. using filter 
subplot(2,1,2)
stem(0:127,impresponse,'red')
grid on 
xlabel('n') 
ylabel('y[n]')
title('Simulated impulse response using filter function')
% first 20 ponint visualization 
figure
subplot(2,1,1) 
stem(0:19,analres(1:20))
grid on 
xlabel('n') 
ylabel('y[n]')
title('Analytical impulse response visualization - first 20 samples')
subplot(2,1,2)
stem(0:19,impresponse(1:20),'red')
grid on 
xlabel('n') 
ylabel('y[n]')
title('Simulated impulse response using filter function')
%% Ex. 1.2
clear 
figure
% a : simulation
n = -10 : 1 : 100 ; % define n 
imp = zeros(length(n)) ; % impulse pre-allocate
imp(n==0) = 1 ;  % our impulse
% stem impulse , x[n] 
subplot(2,1,1) 
stem(n,imp) 
grid on 
xlabel('n') 
ylabel('x[n]') 
title(' x[n] = \delta [n]')
xlim([-10 100])
% calc. h[n] using filter
a = [ 1 -1.8*cos(pi/16) 0.81 ] ; % y coefs
b = [ 1 0.5 ] ; % x coefs
y = filter(b,a,imp) ; % y = h[n] 
subplot(2,1,2) 
stem(n,y) 
grid on 
xlabel('n') 
ylabel('h[n]') 
title('impulse response using filter func')
xlim([-10 100])
% b , analytical answer
analres = zeros(length(n)); % preallocate analytic response
A1 = (1+5/9*exp(-pi/16*1i))/(1 - exp(-pi/8*1i));
A2 = (1+5/9*exp(pi/16*1i))/(1 - exp(pi/8*1i));
analres(n>=0) = A1*(0.9*exp(1i*pi/16)).^n(n>=0)+...
    A2*(0.9*exp(-1i*pi/16)).^n(n>=0); % analytic response definition using 
                                    % relations derived in report file
figure
% simulation result
subplot(2,1,1) 
stem(n,y,'black') 
grid on 
xlabel('n') 
ylabel('h[n]') 
title(' impulse response using filter func')
xlim([-10 100])
% analytica result
subplot(2,1,2) 
stem(n, analres, 'm') 
grid on 
xlabel('n') 
ylabel('h[n]') 
title('Analytical impulse response')
xlim([-10 100])
%% Ex 1.3
clear 
figure
% a 
% equation coefs : 
a = [ 1 -1.8*cos(pi/16) 0.81 ] ; % characteristic polynomial
roots = roots(a) % roots  
% p1n^n u[n]
n = -10 : 100 ; 
p1n=zeros(1,length(n));
p1n(n>=0)= roots(1).^n(n>=0); 
% real of p1^n*u[n]
subplot(2,1,1)
stem(n,real(p1n),'cyan')
xlim([-10 100])
grid on 
xlabel('n')
ylabel('Real Part')
title('real part of  p_{1}^{n} u[n]')
% imag of p1^n*u[n]
subplot(2,1,2)
stem(n,imag(p1n),'green')
grid on 
xlabel('n')
ylabel('Image part')
title('Image part of  p_{1}^{n} u[n]')
xlim([-10 100])
% p2^n u[n]
figure
p2n=zeros(1,length(n));
p2n(n>=0)= roots(2).^n(n>=0); 
% real of p2^n*u[n]
subplot(2,1,1)
stem(n,real(p2n),'cyan')
xlim([-10 100])
grid on 
xlabel('n')
ylabel('Real Part')
title('real part of  p_{2}^{n} u[n]')
% imag of p2^n*u[n]
subplot(2,1,2)
stem(n,imag(p2n),'green')
grid on 
xlim([-10 100])
xlabel('n')
ylabel('Image part')
title('Image part of p_{2}^{n} u[n]')
% b 
% solving equation for alpha and beta
A = [1 1 ; roots(1) roots(2) ] ;
B = [1 ;0.5+1.8*cos(pi/16) ] ;
answers = A\B ; 
alpha = answers(1) 
beta = answers(2) 
h = alpha * p1n + beta * p2n ; % h[n]   imag(h)= 0(epsislon) in matlab
h = real(h);                             % so real(h) = h
figure
subplot(2,1,1)
stem(n,h,'m')
grid on
xlabel('n')
ylabel('h[n]')
title('exercise 1.3 part b result')
xlim([-10 100])
% part 1.2 a results - to compare and verify - 
clear 
n = -10 : 1 : 100 ; % define n 
imp = zeros(length(n)) ; % impulse pre-allocate
imp(n==0) = 1 ;  % our impulse
a = [ 1 -1.8*cos(pi/16) 0.81 ] ; % y coefs
b = [ 1 0.5 ] ; % x coefs
y = filter(b,a,imp) ; % y = h[n] 
subplot(2,1,2)
stem(n,y,'cyan')
grid on
xlabel('n')
ylabel('h[n]')
title('exercise 1.2 part a result - simulation with filter func.')
xlim([-10 100])








