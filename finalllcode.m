%% part 1 : Frequemcy domain (dual tone) :
clear all; clc;
%% Generation dual tone:
Fs = 16000; F1 = 3000; F2 = 5000;
t = (0:Fs-1)/Fs;                                       %range of time
y = sin(2*pi*F1*t) + sin(2*pi*F2*t);                                      
%% Fourier transform & spectrum & spectral density:
N = length(y);
FT = fftshift(fft(y));                                %fourier transforim of signal
f=linspace(-Fs/2,Fs/2,length(FT));                    %range of frequency
%%%%%%% plotting
figure;
subplot(2,1,1); plot(f,abs(FT)/Fs)
title('Spectrum in Frequency domain');
xlabel('Frequency (Hz)');
ylabel('magnitude');
psd=(abs(FT)/Fs).^2;                                  %power spectral density
subplot(2,1,2); plot(f,psd)          
title('spectral density');
xlabel('Frequency (Hz)');
ylabel('dBW/bin');
%% Apply a 3rd order LPF:
Fc = 4000;
[q , r] = butter(3, Fc/(Fs/2));                              % design low pass filter
figure; freqz(q , r);
fsignal = filter(q , r , y);
FT1 = fftshift(fft(fsignal));                                %fourier transforim of filtered signal
figure; plot(f,abs(FT1)/Fs)
title('Spectrum of filtered signal in Frequency domain');
xlabel('Frequency (Hz)');
ylabel('magnitude');
%% using voice signal:
[p,fs]= audioread('eric.wav');         %upload the voice
p=p(1:5*fs);
sound(p,fs)
ts = 1/fs ;
Ns = length(p);
T=fftshift(fft(p));                    %fourier transform of voice signal
f1=linspace(-fs/2,fs/2,length(T));     %range of frequency domain
%%%%%%%%%%%% plotting
figure; 
subplot(2,1,1); plot(f1,abs(T)/fs);
title('spectrum of voice in frequency domain');
xlabel('Frequency (Hz)');
ylabel('magnitude');
psd2=(abs(T)/fs).^2;                  %power spectral density of voice signal
subplot(2,1,2); plot(f1,psd2)
title('spectral density of voice');
xlabel('Frequency (Hz)');
ylabel('dBW/bin');
%%%filtered voice:
[l , k] = fir1(3, Fc/(fs/2));
fsignal1 = filter(l , k , p);
pause(6)
sound(fsignal1,fs)
pause(6)
FT2 = fftshift(fft(fsignal1));                           %fourier transforim of filtered voice signal
figure; plot(f1,abs(FT2)/fs)
title('Spectrum of filtered voice in Frequency domain');
xlabel('Frequency (Hz)');
ylabel('magnitude');

%% part 2 : time domain (triple tone)
clear all; clc
%% generate triple tone
Fs = 20000;
t = (0:Fs-1)/Fs;
F1=100;
F2=2000;
F3=7000;
y = 5*sin(2*pi*F1*t) +5*sin(2*pi*F2*t) +5*sin(2*pi*F3*t);
%%%%%%%%%%%%% plotting
figure;
plot(t, y);
nfft=length(y);
nfft2=2.^nextpow2(nfft);
fy=fft(y,nfft2);                     %fourier transform of signal
fy=fy(1:nfft2/2);
xfft= Fs.*(0:nfft2/2-1)/nfft2;       % range of frequency 
figure;
plot(xfft,abs(fy)/Fs);

%% applying low pass filter with its coefficient Y(t)

coff=[1.0000 0.7303 0.5334 0.3895 0.2845 0.2077 0.1517 0.1108 0.0809 0.0591 0.0432 0.0315 0.023 0.0168 0.0123 0.009 0.0065 0.0048 0.0035 0.0026 0.0019 0.0014 0.001 0.0007];
figure;
freqz(coff,1,512,Fs);                 %plot of filter
x=filter(coff,1,y);
figure; 
plot(t,x);                            %plot of signal in time domain

%% spectrum of filtered signal
nfft9=length(x);
nfft4=2.^nextpow2(nfft9);
fx=fft(x,nfft4);                               %fourier transform of filtered signal
fx=fx(1:nfft4/2);
xfft3= Fs.*(0:nfft4/2-1)/nfft4;                %range of frequency
figure;
plot(xfft3,abs(fx)/Fs);
