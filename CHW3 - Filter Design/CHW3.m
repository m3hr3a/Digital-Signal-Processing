%% HW7 - Simulation Part 
%% DSP - Dr. BabaeiZadeh 
%% Mehrsa Pourya 95101247
%% 7.5 
clear ; 
clc ; 
figure ; 
M = 90 ; 
n = 0 : M ; 
%delta = 0.01 ; 
%A = 40 ; 
beta = 3.395 ; 
%dw = 0.05 * pi ; 
alpha = M/2 ; 
hd = (sin(0.625 * pi *(n - 45)) - sin(0.3 * pi * (n - 45))) ./ (pi * (n - 45)) ; 
hd(n==45) = 0.625 - 0.3 ; 
wk = besseli(0,beta * (1 - ((n-alpha) / alpha).^2).^ (1/2)) ./ besseli(0,beta) ; 
h = hd.*wk ; 
N = (M+1) * 10; 
w = linspace(-pi, pi - 2*pi/N , N) ; 
w = w(ceil((1+length(w))/2):end); 
H = fftshift(fft(h,N)); 
H = H(ceil((1+length(H))/2):end); 
stem(n,h) 
xlabel('n')
ylabel('h[n]')
title('7.5 h[n]')
grid on 
figure
plot(w/pi,abs(H))
title('7.5 , h[n] = w[n] * hd[n] , M = 90') 
xlabel('\omega / pi') 
ylabel('|H(e^{jw})|')
grid minor
hold on
[a,b] = findpeaks(abs(H)) ; 
scatter(w(b(13))/pi,abs(H(b(13))),'fill','red')
text(w(b(12)+7)/pi,abs(H(b(13))+0.05),['SB OvSh=',num2str(a(13))])
scatter(w(b(20))/pi,abs(H(b(20))),'fill','blue')
text(w(b(20)+7)/pi,abs(H(b(20))+0.05),['PB OvSh=',num2str(a(20)-1)])
scatter(w(b(21))/pi,abs(H(b(21))),'fill','green')
text(w(b(21)+5)/pi,abs(H(b(21))+0.05),['SB OvSh=',num2str(a(21))])
hold on
[a,b] = findpeaks(1-abs(H)) ; 
scatter(w(b(15))/pi,abs(H(b(15))),'fill','m')
text(w(b(15)+0)/pi,abs(H(b(15))-0.05),['PB UnSh=',num2str(a(15))])
legend('|H(e^{jw})|','StopBabd Overshoot','PassBabd Overshoot',...
    'StopBabd Overshoot','PassBabd Undershoot','Location','best')
%% 7-6 
clear ; 
clc ; 
figure ; 
delta = 0.03 ; 
A = -20*(log(delta)/log(10)) ;
beta = 0.5842 * (A - 21) ^ (0.4)+0.07886 * (A - 21) ;
dw = 0.05 * pi ; 
M = ceil( (A - 8) / (beta * dw) ) ; 
n = 0 : M ; 
alpha = M/2 ; 
hd = (sin(0.3 * pi *(n - alpha)) - 2*sin(0.5 * pi * (n - alpha))) ./ (pi * (n - alpha)); 
hd(n == 33) = 0.3 - 1 + 2 ; 
wk = besseli(0,beta * (1 - ((n-alpha) / alpha).^2).^ (1/2)) ./ besseli(0,beta) ; 
h = hd.*wk ; 
N = (M+1) * 10; 
w = linspace(-pi, pi - 2*pi/N , N) ; 
H = fftshift(fft(h,N)); 
w = w(ceil((1+length(w))/2):end); 
H = H(ceil((1+length(H))/2):end); 
stem(n,h) 
xlabel('n')
ylabel('h[n]')
title('7.6 h[n]')
xlim([0 M])
grid on 
figure
plot(w/pi,abs(H))
title('7.6 , h[n] = w[n] * hd[n] , M = 66') 
xlabel('\omega / pi') 
ylabel('|H(e^{jw})|')
grid minor
hold on
[a,b] = findpeaks(abs(H)) ; 
scatter(w(b(5))/pi,abs(H(b(5))),'fill','red')
text(w(b(5)+7)/pi,abs(H(b(5))+0.05),['PB1 OvSh=',num2str(a(5)-1)])
scatter(w(b(10))/pi,abs(H(b(10))),'fill','blue')
text(w(b(10)+7)/pi,abs(H(b(10))+0.05),['SB OvSh=',num2str(a(10))])
scatter(w(b(11))/pi,abs(H(b(11))),'fill','green')
text(w(b(11)+5)/pi,abs(H(b(11))+0.05),['PB2 OvSh=',num2str(a(11)-2)])
hold on
[a,b] = findpeaks(1-abs(H)) ; 
scatter(w(b(4))/pi,abs(H(b(4))),'fill','m')
text(w(b(4)+0)/pi,abs(H(b(4))-0.1),['PB1 UnSh=',num2str(a(4))])
[a,b] = findpeaks(2-abs(H)) ; 
scatter(w(b(11))/pi,abs(H(b(11))),'fill','c')
text(w(b(11)+0)/pi,abs(H(b(11))+0.2),['PB2 UnSh=',num2str(a(11))])