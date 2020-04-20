%% DSP - HW10 
%% DTMF 
%% Instructor : Dr. M. BabaeiZadeh
%% Student : Mehrsa Pourya
%% load data 
clear 
clc
close all 
% cd to place of needed files
load('All_Filters.mat')                  % given filters
load('h8.mat')                           % designed filter, passband 1633
load('hlow.mat')                         % lowpass filter 200hz
load('bandpass.mat')                     % bandpass filter i use for noise
                                         % reduction 
[dsnonoise,fs] = audioread('DialedSequence_NoNoise.wav');
[dssnr0,~] = audioread('DialedSequence_SNR00dB.wav');
[dssnr10,~] = audioread('DialedSequence_SNR10dB.wav');
[dssnr20,~] = audioread('DialedSequence_SNR20dB.wav');
[dssnr30,~] = audioread('DialedSequence_SNR30dB.wav');
%% filters'  frequency response
N = 100 * max([length(h1),length(h2),length(h3),length(h4),length(h5),...
    length(h6),length(h7)]);          % N : #points dft
w = 2* pi * (0 : N-1) /N ; 
mid = ceil(length(w)/2) + 1; 
w(mid:end) = w(mid:end) - 2 * pi ; 
f = fftshift(w) * fs / 2 / pi ; 
H(1,:)=abs(fftshift(fft(h1,N)));
H(2,:)=abs(fftshift(fft(h2,N)));
H(3,:)=abs(fftshift(fft(h3,N)));
H(4,:)=abs(fftshift(fft(h4,N)));
H(5,:)=abs(fftshift(fft(h5,N)));
H(6,:)=abs(fftshift(fft(h6,N)));
H(7,:)=abs(fftshift(fft(h7,N)));
H(8,:)=abs(fftshift(fft(h8,N)));
for i = 1 : 8
    plot(f,H(i,:))
    [a,b] = max(H(i,:));
    hold on
    scatter(f(b),a,'fill')
    hold on
    text(f(b),a,['h',num2str(i)])
    text(-400,1-0.1*i,['f(h',num2str(i),') = ',num2str(f(b))])
end
xlim([-1800 1800])
grid minor
xlabel('f(Hz)')
ylabel('|H(e^{jw})|')
title('Magnitude of frequency response of filters')
clear a b H i mid w 
figure 
plot(f,abs(fftshift(fft(hlow,N))))
xlim([-800 800])
grid minor
xlabel('f(Hz)')
ylabel('|H(e^{jw})|')
title('Magnitude of frequency response of lowpass(200Hz-500Hz)')
%% no noise dialed sequence decode
figure
signal = dsnonoise;
z = zeros(10,1);
t = (0 : length(signal)-1)/fs ; 
plot(t,signal)
hold on
% first step is to seperate each dial section , i find consecutive zeros 
% and then find start and end of each section according to zero intervals
zplaces = find(signal == 0);                        
ed = zplaces - [0 ; zplaces(1:end-1)]; 
ed = find(ed ~= 1);                       % non consecutive zeros starter
k = 1 ; 
% a zero starter is a section's end if many samples after it is zero :
% de(i) : zero interval(i)'s start = dial section(i)'s end
for i = 1 : length(ed)
    if(isequal(signal(zplaces(ed(i)):zplaces(ed(i))+9),z))
        dend(k) = t(zplaces(ed(i))-1);
        de(k) = zplaces(ed(i))-1;
        k = k+1 ;
    end
end
scatter(dend,zeros(length(dend),1),'fill','g')
st = [zplaces(2:end) ; 0] - zplaces ; 
st = find(st ~= 1);                      % non consecutive zeros ender            
k = 1 ; 
% a zero ender is a section's start if many samples before it is zero :
% ds(i) : zero interval(i)'s end = dial section(i)'s start
for i = 1 : length(st)-1
    if(isequal(signal(zplaces(st(i)-9):zplaces(st(i))),z))
       dst(k) = t(zplaces(st(i))+1);
       ds(k) = zplaces(st(i))+1;
       k = k+1 ;
    end
end
scatter(dst,zeros(length(dst),1),'fill','r')
title('dialed sequence no noise ,for each section, green : end , red : start')
xlabel('t')
xlim([0 t(end)])
clear zplaces z signal k i st ed dst dend
clc
% decoder of each section , find row using low freq ans col using high freq
series = []; 
for i = 1 : length(ds)
sec = dsnonoise(ds(i):de(i));
a1l = conv(conv(sec,h1).^2,hlow);
a2l = conv(conv(sec,h2).^2,hlow);
a3l = conv(conv(sec,h3).^2,hlow);
a4l = conv(conv(sec,h4).^2,hlow);
[~,row] = max([mean(a1l) mean(a2l) mean(a3l) mean(a4l)]);
a1h = conv(conv(sec,h5).^2,hlow);
a2h = conv(conv(sec,h6).^2,hlow);
a3h = conv(conv(sec,h7).^2,hlow);
a4h = conv(conv(sec,h8).^2,hlow);
[~,col] = max([mean(a1h) mean(a2h) mean(a3h) mean(a4h)]);
series =[series,' ',whichdigit(row,col)];
end
ds_no_noise_decoded = series
clear a1l a2l a3l a4l a1h a2h a3h a4h series sec i 
%% snr 30
% for noisy signals seperating sections need another algorithm 
figure
s = conv(dssnr30,bandpass);             % we filter signal to reduce noise
zsample = ceil(0.2*fs);                 % number of zero interval samples 
t = (0 : length(s)-1) / fs ; 
% i convolve asb(signal) with a 1 pulse(length=zero samples) and local
% minimums of result corrospond to end of each dial section and +zero
% samples -> start section , i repeat this algorithm for all noisy signals
u = ones(zsample,1);
p = abs(conv(u,abs(s)));
[a b] = findpeaks(1-p,'MinPeakDistance',400);
plot(p)
hold on
bn = find(a>-16);
bn(1)=[];
scatter(b(bn),p(b(bn)),'fill')
title('dialed sequence SNR = 30, convolution with ones(1,0.2*fs)')
grid on
xlabel('samples')
figure
plot(t,s)
hold on 
title('dialed sequence SNR = 30 ,for each section, green : end , red : start')
xlabel('t')
xlim([0 t(end)])
de = b(bn)-zsample;
ds = [zsample;b(bn(1:end-1))];
scatter(t(ds),zeros(1,length(ds)),'fill','r')
scatter(t(de),zeros(size(de)),'fill','g')
series = []; 
for i = 1 : length(ds)
sec = s(ds(i):de(i));
a1l = conv(conv(sec,h1).^2,hlow);
a2l = conv(conv(sec,h2).^2,hlow);
a3l = conv(conv(sec,h3).^2,hlow);
a4l = conv(conv(sec,h4).^2,hlow);
[~,row] = max([mean(a1l) mean(a2l) mean(a3l) mean(a4l)]);
a1h = conv(conv(sec,h5).^2,hlow);
a2h = conv(conv(sec,h6).^2,hlow);
a3h = conv(conv(sec,h7).^2,hlow);
a4h = conv(conv(sec,h8).^2,hlow);
[~,col] = max([mean(a1h) mean(a2h) mean(a3h) mean(a4h)]);
series =[series,' ',whichdigit(row,col)];
end
ds_snr30_decoded = series
clear a1l a2l a3l a4l a1h a2h a3h a4h series sec i 
%% snr 20 
figure
s =conv(dssnr20,bandpass);
zsample = ceil(0.2*fs);
t = (0 : length(s)-1) / fs ; 
sstop = ceil(0.2*fs) ; 
u = ones(sstop,1);
p = abs(conv(u,abs(s)));
[a b] = findpeaks(1-p,'MinPeakDistance',400);
plot(p)
hold on
bn = find(a>-30);
bn(1)=[];
scatter(b(bn),p(b(bn)),'fill')
title('dialed sequence SNR = 20, convolution with ones(1,0.2*fs)')
grid on
xlabel('samples')
figure
plot(t,s)
title('dialed sequence SNR =20 ,for each section, green : end , red : start')
xlabel('t')
xlim([0 t(end)])
hold on
de = b(bn)-zsample;
ds = [zsample+450;b(bn(1:end-1))];
scatter(t(ds),zeros(1,length(ds)),'fill','r')
scatter(t(de),zeros(size(de)),'fill','g')
series = []; 
for i = 1 : length(ds)
sec = s(ds(i):de(i));
a1l = conv(conv(sec,h1).^2,hlow);
a2l = conv(conv(sec,h2).^2,hlow);
a3l = conv(conv(sec,h3).^2,hlow);
a4l = conv(conv(sec,h4).^2,hlow);
[~,row] = max([mean(a1l) mean(a2l) mean(a3l) mean(a4l)]);
a1h = conv(conv(sec,h5).^2,hlow);
a2h = conv(conv(sec,h6).^2,hlow);
a3h = conv(conv(sec,h7).^2,hlow);
a4h = conv(conv(sec,h8).^2,hlow);
[~,col] = max([mean(a1h) mean(a2h) mean(a3h) mean(a4h)]);
series =[series,' ',whichdigit(row,col)];
end
ds_snr20_decoded = series
clear a1l a2l a3l a4l a1h a2h a3h a4h series sec i 
%% snr 10
figure
s =conv(dssnr10,bandpass);
zsample = ceil(0.2*fs);
t = (0 : length(s)-1) / fs ; 
sstop = ceil(0.2*fs) ; 
u = ones(sstop,1);
p = abs(conv(u,abs(s)));
[a b] = findpeaks(1-p,'MinPeakDistance',400);
plot(p)
hold on
bn = find(a>-65);
bn(1)=[];
scatter(b(bn),p(b(bn)),'fill')
title('dialed sequence SNR = 10, convolution with ones(1,0.2*fs)')
xlabel('t')
grid on
xlabel('samples')
figure
plot(t,s)
title('dialed sequence SNR =10 ,for each section, green : end , red : start')
xlabel('t')
hold on
de = b(bn)-zsample;
ds = [zsample+450;b(bn(1:end-1))];
scatter(t(ds),zeros(1,length(ds)),'fill','r')
scatter(t(de),zeros(size(de)),'fill','g')
series = []; 
for i = 1 : length(ds)
sec = s(ds(i):de(i));
a1l = conv(conv(sec,h1).^2,hlow);
a2l = conv(conv(sec,h2).^2,hlow);
a3l = conv(conv(sec,h3).^2,hlow);
a4l = conv(conv(sec,h4).^2,hlow);
[~,row] = max([mean(a1l) mean(a2l) mean(a3l) mean(a4l)]);
a1h = conv(conv(sec,h5).^2,hlow);
a2h = conv(conv(sec,h6).^2,hlow);
a3h = conv(conv(sec,h7).^2,hlow);
a4h = conv(conv(sec,h8).^2,hlow);
[~,col] = max([mean(a1h) mean(a2h) mean(a3h) mean(a4h)]);
series =[series,' ',whichdigit(row,col)];
end
ds_snr10_decoded = series
clear a1l a2l a3l a4l a1h a2h a3h a4h series sec i 
%% snr 0
figure
s =conv(dssnr0,bandpass);
zsample = ceil(0.2*fs);
t = (0 : length(s)-1) / fs ; 
sstop = ceil(0.2*fs) ; 
u = ones(sstop,1);
p = abs(conv(u,abs(s)));
[a b] = findpeaks(1-p,'MinPeakDistance',400);
plot(p)
hold on
bn = find(a>-128);
bn(1)=[];
scatter(b(bn),p(b(bn)),'fill')
title('dialed sequence SNR = 0, convolution with ones(1,0.2*fs)')
grid on
xlabel('samples')
figure
plot(t,s)
title('dialed sequence SNR =0 ,for each section, green : end , red : start')
xlabel('t')
hold on
de = b(bn)-zsample;
ds = [zsample+450;b(bn(1:end-1))];
scatter(t(ds),zeros(1,length(ds)),'fill','r')
scatter(t(de),zeros(size(de)),'fill','g')
series = []; 
for i = 1 : length(ds)
sec = s(ds(i):de(i));
a1l = conv(conv(sec,h1).^2,hlow);
a2l = conv(conv(sec,h2).^2,hlow);
a3l = conv(conv(sec,h3).^2,hlow);
a4l = conv(conv(sec,h4).^2,hlow);
[~,row] = max([mean(a1l) mean(a2l) mean(a3l) mean(a4l)]);
a1h = conv(conv(sec,h5).^2,hlow);
a2h = conv(conv(sec,h6).^2,hlow);
a3h = conv(conv(sec,h7).^2,hlow);
a4h = conv(conv(sec,h8).^2,hlow);
[~,col] = max([mean(a1h) mean(a2h) mean(a3h) mean(a4h)]);
series =[series,' ',whichdigit(row,col)];
end
ds_snr0_decoded = series
clear a1l a2l a3l a4l a1h a2h a3h a4h series sec i 
%% tabel
function digit = whichdigit(row,col)
    temp = (col-1) + 4 * (row-1) ; 
    switch(temp)
        case 0
            digit = '1';
        case 1
            digit = '2';
        case 2
            digit = '3';
        case 3
            digit = 'A';
        case 4
            digit = '4';
        case 5
            digit = '5';
        case 6
            digit = '6';
        case 7
            digit = 'B';
        case 8
            digit = '7';
        case 9
            digit = '8';
        case 10
            digit = '9';
        case 11
            digit = 'C';
        case 12
            digit = '*';
        case 13
            digit = '0';
        case 14
            digit = '#';
        case 15
            digit = 'D';
    end
end