%% Discrete-Time Signal Processing 
%% Dr Babaei-Zade Spring 2019 
%% CHW2
%% Student : Mehrsa Pourya 95101247
%% Project 1 
%% Ex 1.1 . b 
clc;clear;close all;
L = 12 ; 
n = 0 : L-1 ; 
r(1:L) = 1 ; 
for m = 5 : 10 
    figure
    subplot(3,2,[1 2]) 
    stem(n,r,'.','Linewidth',2)
    grid minor
    xlabel('n')
    ylabel('r[n]')
    title(['Ex 1.1 b -- r[n] in time domain (L=12),then we calculate dtft,N=L*',num2str(m)])
    xlim([0 L-1])
    N = m * L ; % m = N/L
    [R,w] = dtft(r,0,N); 
    subplot(3,2,3) 
    plot(w,real(R),'m','Linewidth',0.8)
    grid minor
    xlabel('w')
    ylabel('Re[R(e^({jw})]')
    title('real part of dtft of r[n]')
    xlim([-pi pi])
    subplot(3,2,4) 
    plot(w,imag(R),'g','Linewidth',0.8)
    grid minor
    xlabel('w')
    ylabel('Im[R(e^({jw})]')
    title('Image part of dtft of r[n]')
    xlim([-pi pi])
    subplot(3,2,[5 6]) 
    plot(w,abs(R),'m','Linewidth',0.8)
    grid minor
    xlabel('w')
    ylabel('abs[R(e^({jw})]')
    title('magnitude of dtft of r[n]')
    xlim([-pi pi])
end
%% Ex 1.1 . c
clear;
for L = [12 15] 
    n = 0 : L-1 ; 
    r(1:L) = 1 ; 
    for m = 8
        figure
        subplot(6,1,1) 
        stem(n,r,'.','Linewidth',2)
        grid minor
        xlabel('n')
        ylabel('r[n]')
        title(['Ex 1.1 c -- r[n] in time domain (L=15),then we calculate dtft , N=L*',num2str(m)])
        xlim([0 L-1])
        N = m * L ; % m = N/L
        [R,w] = dtft(r,0,N); 
        subplot(6,1,[2 3 4])
        zeroplaces = find(abs(R) < 1e-4);
        wzero = w(zeroplaces);
        zerowdist = wzero(2:end) - wzero(1:end-1);
        plot(w,abs(R),'m','Linewidth',0.8)
        hold on
        scatter(wzero,zeros(size(wzero)),'fill')
        hold on 
        scatter(w(abs(R)==max(abs(R))),max(abs(R)),'fill')
        text(w(abs(R)==max(abs(R)))-1,max(abs(R))+2,['peak value =',num2str(max(R))])
        for z = 1 : length(wzero)
            text(wzero(z),R(w==wzero(z))-2,['#',num2str(z)])
        end
        grid minor
        xlabel('w')
        ylabel('|R(e^({jw})|')
        legend('|R(e^({jw})|','zero number(#)','peak/dc','Location','Best')
        title('magnitude of dtft of r[n]')
        ylim([-4 max(R)+4])
        xlim([-pi pi])
        subplot(6,1,[5 6])
        stem(2:length(wzero),zerowdist,'.black','Linewidth',0.8)
        grid minor
        xlabel('i')
        title(['R(e^{jw}) has ',num2str(length(wzero)),...
            ' zeros and the distanse between ith and (i-1)th zero is plotted'])
        xlim([2 length(wzero)])
    end
end
%% Ex 1.2
clear;figure;
L = 12 ; 
N = 8 * 12; 
n = 0 : L - 1 ; 
r = ones(1,L); 
[R1,w] = dtft(r,0,N); 
R2 = asinc(L,w) .* exp(-1i*w*(L-1)/2) ; 
subplot(2,2,[1 2])
stem(n,r,'fill')
xlim([0 L-1])
title('Ex 1.2 -- r(n), L = 12 , N = 5 * 12 = 60')
xlabel('n')
grid minor
subplot(2,2,3)
plot(w,abs(R1))
title('|R(e^{jw})| using dtft function')
xlabel('w')
xlim([-pi pi])
grid minor
subplot(2,2,4)
plot(w,abs(R2))
title('|R(e^{jw})| using asinc and theory')
xlabel('w')
xlim([-pi pi])
grid minor
%% Ex 1.3 a b
clc;figure;
L = 12 ; 
N = L * 50 ; 
n = 0 : L-1 ; 
r(1:L) = 1 ; 
subplot(7,1,1) 
stem(n,r','.','Linewidth',2)
grid minor
xlabel('n')
ylabel('r[n]')
title(['Ex 1.3 a,b -- r[n] in time domain (L=12),then we calculate dtft , N=',num2str(N)])
xlim([0 L-1])
[R,w] = dtft(r,0,N); 
subplot(7,1,[2 3]) 
plot(w,abs(R),'m','Linewidth',0.9)
grid minor
xlabel('w')
ylabel('abs[R(e^({jw})]')
title('magnitude of dtft of r[n]')
xlim([-pi pi])
subplot(7,1,[4 5]) 
plot(w,angle(R),'g','Linewidth',0.9)
grid minor
xlabel('w')
ylabel('phase[R(e^({jw})]')
title('phase of dtft of r[n] using angle')
xlim([-pi pi])
subplot(7,1,[6 7]) 
plot(w,unwrap(angle(R)),'g','Linewidth',0.9)
grid minor
xlabel('w')
ylabel('phase[R(e^({jw})]')
title('unwrapped phase of dtft of r[n] using angle')
xlim([-pi pi])
%% Project 2 Ex 2.3 a & b
% a 
clear;figure;
n = 0 : 20 ; 
x = 0.9.^n .* cos(2*pi*n/sqrt(31));
[X,w] = dtft(x,0,21*5) ; 
subplot(3,2,1)
plot(w,abs(X))
grid minor
ylabel('|X(e^{jw})|')
xlabel('w')
title('Ex. 2.3 a  -- plot 1(p1)')
xlim([-pi pi])
subplot(3,2,2)
plot(w,unwrap(angle(X)))
grid minor
ylabel('<[X(e^{jw})]')
xlabel('w')
xlim([-pi pi])
title('p2')
[Xflip,wf]= flipDTFT(X,w);
subplot(3,2,3)
plot(wf,abs(Xflip))
grid minor
ylabel('|X(e^{-jw})|')
xlabel('w')
xlim([-pi pi])
title('p3')
subplot(3,2,4)
plot(wf,unwrap(angle(Xflip)))
grid minor
ylabel('<[X(e^{-jw})]')
xlabel('w')
xlim([-pi pi])
Xconj = conj(X);
title('p4')
subplot(3,2,5)
plot(w,abs(Xconj))
grid minor
ylabel('|X*(e^{jw})|')
xlabel('w')
xlim([-pi pi])
title('p5')
subplot(3,2,6)
plot(w,unwrap(angle(Xconj)))
grid minor
ylabel('<[X*(e^{jw})]')
xlabel('w')
xlim([-pi pi])
title('p6, p6=p4 & p3=p5 --> CS of DTFT[real signal]')
% b
figure
y = 1i * x ; 
[Y,w] = dtft(y,0,21*5) ; 
subplot(4,2,1)
plot(w,abs(Y))
grid minor
ylabel('|Y(e^{jw})|')
xlabel('w')
title('Ex. 2.3 b  -- plot 1(p1)')
xlim([-pi pi])
subplot(4,2,2)
plot(w,unwrap(angle(Y)))
grid minor
ylabel('<[Y(e^{jw})]')
xlabel('w')
xlim([-pi pi])
title('p2')
[Yflip,wf]= flipDTFT(Y,w);
subplot(4,2,3)
plot(wf,abs(Yflip))
grid minor
ylabel('|Y(e^{-jw})|')
xlabel('w')
xlim([-pi pi])
title('p3')
subplot(4,2,4)
plot(wf,unwrap(angle(Yflip)))
grid minor
ylabel('<[Y(e^{-jw})]')
xlabel('w')
xlim([-pi pi])
Yconj = conj(Y);
title('p4')
subplot(4,2,5)
plot(w,abs(Yconj))
grid minor
ylabel('|Y*(e^{jw})|')
xlabel('w')
xlim([-pi pi])
title('p5')
subplot(4,2,6)
plot(w,unwrap(angle(Yconj)))
grid minor
ylabel('<[Y*(e^{jw})]')
xlabel('w')
xlim([-pi pi])
title('p6')
subplot(4,2,[7 8])
plot(w,unwrap(angle(Yflip))-unwrap(angle(Yconj)))
title('p7,p4-p6=\pi & p3=p5 -> Y(e^{-jw})=e^{j\pi}Y*(e^{jw})=-Y*(e^{jw}) -> CAS of DTFT(IMG sig.)')
xlim([-pi pi])
ylim([-1 pi+1])
grid minor
xlabel(w)
ylabel('p4-p6')
%% c 
clear;figure; 
n = -29 : 29 ; 
ve = exp(1i*2*pi*n.^2/25); 
subplot(2,2,1)
stem(n,real(ve),'.')
grid minor
ylabel('real[ve(n)]')
xlabel('n')
title('Ex. 2.3 c  -- plot 1(p1)')
xlim([-29 29])
subplot(2,2,2)
stem(n,imag(ve),'.')
grid minor
ylabel('imag[ve(n)]')
xlabel('n')
xlim([-29 29])
title('p2 , p1 & p2 : both even -> ve[n] is even')
[Ve,w] = dtft(ve,29, 59*5) ; 
subplot(2,2,3)
plot(w,real(Ve))
grid minor
ylabel('real[Ve(e^{jw})]')
xlabel('n')
title('p3')
xlim([-pi pi])
subplot(2,2,4)
plot(w,imag(Ve))
grid minor
ylabel('imag[Ve(e^{jw})]')
xlabel('w')
xlim([-pi pi])
title('p4, p3 & p4 :both even -> DTFT(ve[n]) is even')
%% d
clear;figure; 
n = -19 : 19 ; 
vo = n; 
subplot(2,2,1)
stem(n,real(vo),'.')
grid minor
ylabel('real[vo(n)]')
xlabel('n')
title('Ex. 2.3 d  -- plot 1(p1)')
xlim([-19 19])
subplot(2,2,2)
stem(n,imag(vo),'.')
grid minor
ylabel('imag[vo(n)]')
xlabel('n')
xlim([-19 19])
title('p2=0 , p1:Odd -> vo[n] is real and Odd')
[Vo,w] = dtft(vo,19, 390) ; 
subplot(2,2,3)
plot(w,real(Vo))
grid minor
ylabel('real[Vo(e^{jw})]')
xlabel('n')
title('p3 = 0 ')
xlim([-pi pi])
subplot(2,2,4)
plot(w,imag(Vo))
grid minor
ylabel('imag[Vo(e^{jw})]')
xlabel('w')
xlim([-pi pi])
title('p4 : Odd -> DTFT(vo[n]): Purly imag. and Odd')
%% e 
clear;figure; 
n = -19 : 19 ; 
go = 1i * n; 
subplot(3,2,1)
stem(n,real(go),'.')
grid minor
ylabel('real[Go(n)]')
xlabel('n')
title('Ex. 2.3 e  -- plot 1(p1)')
xlim([-19 19])
subplot(3,2,2)
stem(n,imag(go),'.')
grid minor
ylabel('imag[go(n)]')
xlabel('n')
xlim([-19 19])
title('p2 , p2:Odd -> go[n] : imag. & Odd')
[Go,w] = dtft(go,19, 390) ; 
subplot(3,2,3)
plot(w,real(Go))
grid minor
ylabel('real[Go(e^{jw})]')
xlabel('n')
title('p3&p4=0 -> DTFT[go[n]] : Odd and real')
xlim([-pi pi])
subplot(3,2,4)
plot(w,imag(Go))
grid minor
ylabel('imag[Go(e^{jw})]')
xlabel('w')
xlim([-pi pi])
title('p4=0')
[ Gf , wf ] = flipDTFT( Go, w );
subplot(3,2,5)
plot(wf,real(Gf))
grid minor
ylabel('real[Go(e^{-jw})]')
xlabel('n')
title('p5')
xlim([-pi pi])
subplot(3,2,6)
plot(wf,imag(Gf))
grid minor
ylabel('imag[Go(e^{-jw})]')
xlabel('w')
xlim([-pi pi])
title('p6=0 , p5=-p3 -> Go(e^{jw})=-Go(e^{-jw})')
%% Ex 3.1 
clear;figure;
n = 0 : 30 ; 
x = 0.9 .^ n ; %system : 1/1- 0.9 e^(-jw)
b = 1 ; 
a = [1 -0.9] ; 
h = zeros(size(n)); 
h(1) = 1 ;   % h is our impulse 
N = 300 ; 
[ X, w ] = freqz(b,a,N,'whole');
mid = ceil(N/2) + 1;
w(mid:N) = w(mid:N) - 2*pi; 
w = fftshift(w);
X = fftshift(X);
subplot(3,2,[1 2])
stem(n,x)
xlim([0 30])
grid on
xlabel('n')
ylabel('x[n]')
title('Ex 3.1 ,p1 x[n] = 0.9^n u[n]')
subplot(3,2,3)
plot(w,abs(X))
xlim([-pi pi])
grid on
xlabel('w')
ylabel('|X(e^{jw})|')
title('p2 mag. ->,Even using freqz')
subplot(3,2,4)
plot(w,unwrap(angle(X)))
xlim([-pi pi])
grid on
xlabel('w')
ylabel('<[X(e^{jw})] using freqz')
title('p3 phase ->,Odd using freqz')
% formula Xtheory = 1/(1-0.9e^(-jw))
Xt = 1./(1-0.9*exp(-1i*w)); 
subplot(3,2,5)
plot(w,abs(Xt))
xlim([-pi pi])
grid on
xlabel('w')
ylabel('|X(e^{jw})|')
title('p4 mag. -> Even,using theory = p2')
subplot(3,2,6)
plot(w,unwrap(angle(Xt)))
xlim([-pi pi])
grid on
xlabel('w')
ylabel('<[X(e^{jw})] using freqz')
title('p5 phase ->  Odd,using theory = p3 ')
%% Ex 4.2 
%% a
clear;figure;
L = 32;
n = 0 : L-1; 
r = ones(size(n)); 
theta = 2 * pi / sqrt(31) ; 
x = r .* exp(1i * theta * n) ; 
[X,w] = dtft(x,0,8*L);
subplot(4,1,1)
stem(n,abs(x),'.') 
grid minor
xlabel('n')
ylabel('|x[n]|')
title('Ex. 4.2 a -- x[n]=r.*e^{j*\theta * n},p1')
xlim([0 L-1])
subplot(4,1,2)
stem(n,unwrap(angle(x)),'.') 
grid minor
xlabel('n')
ylabel('<[x[n]]')
title('p2')
xlim([0 L-1])
subplot(4,1,3)
plot(w,abs(X)) 
hold on 
scatter(w((abs(X)==max(abs(X)))),max(abs(X)),'fill')
text(w((abs(X)==max(abs(X))))+0.1,max(abs(X)),['peak : X = ',...
    num2str(w((abs(X)==max(abs(X))))),' Y = ',num2str(max(abs(X)))])
grid minor
xlabel('w')
ylabel('|X(e^{jw})|')
title('p3 , peak is at 1.129 = theta = 2*pi/sqrt(31)')
xlim([-pi pi])
ylim([ -1 max(abs(X))+5])
subplot(4,1,4)
plot(w,unwrap(angle(X))) 
grid minor
xlabel('n')
ylabel('<[X(e^{jw})]')
title('p4')
xlim([-pi pi])
%% b
clear;figure;
L = 32;
n = 0 : L-1; 
w = 1/2 - 1/2 * cos(2*pi/L*n); 
theta = 2 * pi / sqrt(31) ; 
g = w .* exp(1i * theta * n) ; 
[G,w] = dtft(g,0,8*L);
subplot(4,1,1)
stem(n,abs(g),'.') 
grid minor
xlabel('n')
ylabel('|g[n]|')
title('Ex. 4.2 b -- w[n]=1/2-1/2cos(2\pi n/L), g[n]=w[n].*e^{j*\theta * n},p1')
xlim([0 L-1])
subplot(4,1,2)
stem(n,unwrap(angle(g)),'.') 
grid minor
xlabel('n')
ylabel('<[g[n]]')
title('p2')
xlim([0 L-1])
subplot(4,1,3)
plot(w,abs(G)) 
hold on 
scatter(w((abs(G)==max(abs(G)))),max(abs(G)),'fill')
text(w((abs(G)==max(abs(G))))+0.1,max(abs(G)),['peak : G = ',...
    num2str(w((abs(G)==max(abs(G))))),' Y = ',num2str(max(abs(G)))])
grid minor
xlabel('w')
ylabel('|G(e^{jw})|')
title('p3 , peak is at 1.129 = theta = 2*pi/sqrt(31)')
xlim([-pi pi])
ylim([ -1 max(abs(G))+5])
subplot(4,1,4)
plot(w,unwrap(angle(G))) 
grid minor
xlabel('n')
ylabel('<[G(e^{jw})]')
title('p4')
xlim([-pi pi])
%% Ex 5.1
% a T = 1ms -> fs = 1khz -> max f = 500 hz
% b
clear;figure;
w0 = 2*pi/5 ; 
w = linspace(-pi,pi,1000); 
H = (1-exp(-1i*(w-w0))).*(1-exp(-1i*(w+w0)))./...
    (1-0.9*exp(-1i*(w-w0)))./(1-0.9*exp(-1i*(w+w0)));
subplot(2,1,1) 
plot(w,abs(H))
xlabel('w')
ylabel('|H(e^{jw})|')
grid minor
title('Ex 5.1 a -- p1 , mag.')
xlim([-pi pi])
subplot(2,1,2)
plot(w,(angle(H)))
xlabel('w')
ylabel('<[H(e^{jw})]')
grid minor
xlim([-pi pi])
title('p2 , phaee')
% c w = omega*T =2*pi*60 *(1/1000)=3*pi/25
% d
figure
w0 = 3 * pi / 25 ; 
H = (1-exp(-1i*(w-w0))).*(1-exp(-1i*(w+w0)))./...
    (1-0.9*exp(-1i*(w-w0)))./(1-0.9*exp(-1i*(w+w0)));
bH = [1 -2*cos(w0) 1] ; 
aH = [1 -2*0.9*cos(w0) 0.81] ; 
plot(w,abs(H))
xlabel('w')
ylabel('|H(e^{jw})|')
grid minor
title('Ex 5.1 d -- mag(H) for w0 = 3\pi / 25')
xlim([-pi pi])
% e 
figure
fs = 1000 ; 
t = 0 : 1/fs : 150/fs-1/fs;
mysin = sin(2*pi*60*t) ; 
out = filter(bH,aH,mysin) ; 
subplot(4,1,1)
plot(t,mysin)
xlabel('t')
title('Ex 5.1 e & f , p1 :input =  sin(2\pi*60t)')
grid minor
subplot(4,1,2)
I =fft(mysin)/150;
f = linspace(-fs/2,fs/2,150) ; 
plot(f,fftshift(abs(I)))
xlabel('f')
ylabel('|Input(j\Omega)|')
title('p2 |DTFT(input)|')
ylim([0 0.6])
grid minor
subplot(4,1,3)
plot(t,out)
xlabel('t')
title('p3, output of given filter with wo = 3\pi/2')
hold on
scatter(t(40) , abs(out(40)),'fill')
text(t(41) , abs(out(40))+0.2,['transient time is ',num2str(t(40)),'ms'])
grid minor
subplot(4,1,4)
O= fft(out);
plot(f,fftshift(abs(O))/length(O))
ylim([0 0.6])
xlabel('f')
ylabel('|Output(j\Omega)|')
title('p4 |DTFT(output)|')
grid minor
%% page 26 Ex 1.1
% b 
% mygrpdelay function is attached
% c
% mygrpdlywnstart function is attached
