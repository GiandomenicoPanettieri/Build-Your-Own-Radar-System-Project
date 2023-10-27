clear all;
close all;
clc;
%Constants
fc = 5.8e9;
Tp = 0.1;
N = 10;
c  = 299792458;                   
SamplingFrequency=44100;     


I = audioread("RealPartCWIF.wav")';
I = I-mean(I);
Q = audioread("ImaginaryPartCWIF.wav")';
Q = Q-mean(Q);

AudioDuration = length(I)/SamplingFrequency;  %time duration of input data
% time axis
dt = 1/SamplingFrequency;
t = 0 : dt : AudioDuration-dt;

t_0=2;
t = t(t_0*SamplingFrequency:end);
I = I(t_0*SamplingFrequency:end);
Q = Q(t_0*SamplingFrequency:end);

X = I + 1i * Q;  % Create complex signal


figure(2)
plot(t,angle(X)); xlabel("time [s]"); ylabel("magnitude"); title("CW LOW IF, PHASE(I+jQ) SIGNAL VS TIME PLOT")
%Range Plot
figure(3)
Range = (c/(4*pi*fc))*unwrap(angle(X));
plot(t,-Range); xlabel("time [s]"); ylabel("range[m]"); title("Breathing Pattern CW Low IF ")

%FFT plot
BW = SamplingFrequency;
bin = BW/(length(t)*(N+1));
f = 0:bin:BW-bin;
figure(4)
ZeroPaddingRange = zeros(1,N*length(t));
Range = [Range ZeroPaddingRange];
Range_FFT = fft([Range])/max(abs(fft(Range)));
plot(f,20*log10(abs(Range_FFT))); xlabel("frequency [Hz]"); ylabel("Magnitude"); title("FFT Breathing Pattern");
axis([0 2 -55 0]);
