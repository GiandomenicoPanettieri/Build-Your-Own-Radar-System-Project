clear all;
close all;
clc;
%Constants
fc = 5.8e9;
Tp = 0.1;
N = 4;
c  = 299792458;                   
SamplingFrequency=44100;     

I = audioread("RealPartCWIF.wav")';
I = I-mean(I);
AudioDuration = length(I)/SamplingFrequency;  %time duration of input data
% time axis
dt = 1/SamplingFrequency;
t = 0 : dt : AudioDuration-dt;
%frequency and velocity axis
BandWidth = SamplingFrequency;
bin = SamplingFrequency/(length(I)*(N+1)); 
f = 0 : bin : BandWidth-bin;                        
v = f * (c/(2*fc)); %velocity axis                                 
        
% Ensemble Matrix initialisation
NumberOfTimeWindows = floor(AudioDuration/(Tp));
SamplesInWindow = Tp*SamplingFrequency;         
WindowCollectionMatrixI = zeros(NumberOfTimeWindows, SamplesInWindow*(1+N));

for kk = 1 : NumberOfTimeWindows   
    TemporaryVector = I( ((kk-1)*SamplesInWindow +1) : ((kk)*SamplesInWindow) );
    TemporaryZeroVector = zeros(1, N*length(TemporaryVector));
    TemporaryVectorZeroPadding = [TemporaryVector TemporaryZeroVector];
    WindowCollectionMatrixI(kk,:) = TemporaryVectorZeroPadding; 
end
% FFT 
SpectrogramI = zeros(size(WindowCollectionMatrixI));
for kk = 1 : NumberOfTimeWindows
    SpectrogramI(kk,:) = fft(WindowCollectionMatrixI(kk,:));
end
 SpectrogramI = 20*log10( abs(SpectrogramI) );
 SpectrogramI = CutLowValue(SpectrogramI,-60);
% Normalization row by row
SpectrogramI = SpectrogramI-max(SpectrogramI')';
Q = audioread("ImaginaryPartCWIF.wav")';
Q = Q-mean(Q);
% Ensemble Matrix initialisation
NumberOfTimeWindows = floor(AudioDuration/(Tp));         
WindowCollectionMatrixQ = zeros(NumberOfTimeWindows, SamplesInWindow*(1+N));

for kk = 1 : NumberOfTimeWindows   
    TemporaryVector = Q( ((kk-1)*SamplesInWindow +1) : ((kk)*SamplesInWindow) );
    TemporaryZeroVector = zeros(1, N*length(TemporaryVector));
    TemporaryVectorZeroPadding = [TemporaryVector TemporaryZeroVector];
    WindowCollectionMatrixQ(kk,:) = TemporaryVectorZeroPadding; 
end
% FFT 
SpectrogramQ = zeros(size(WindowCollectionMatrixQ));
for kk = 1 : NumberOfTimeWindows
    SpectrogramQ(kk,:) = fft(WindowCollectionMatrixQ(kk,:));
end
SpectrogramQ = 20*log10( abs(SpectrogramQ) );
SpectrogramQ = SpectrogramQ-max(SpectrogramQ')';
SpectrogramQ = CutLowValue(SpectrogramQ,-10);

X = I + 1i * Q;  % Create complex signal

% Ensemble Matrix initialization for complex signal
WindowCollectionMatrixX = zeros(NumberOfTimeWindows, SamplesInWindow * (1 + N));

% Filling the Ensemble Matrix with Zero Padding for complex signal
for kk = 1 : NumberOfTimeWindows   
    TemporaryVector = X(((kk - 1) * SamplesInWindow + 1) : (kk * SamplesInWindow));
    TemporaryZeroVector = zeros(1, SamplesInWindow * N);
    TemporaryVectorZeroPadding = [TemporaryVector, TemporaryZeroVector];
    WindowCollectionMatrixX(kk, :) = TemporaryVectorZeroPadding; 
end

% FFT for complex signal
SpectrogramX1 = zeros(size(WindowCollectionMatrixX));
for kk = 1 : NumberOfTimeWindows
    SpectrogramX1(kk, :) = fft(WindowCollectionMatrixX(kk, :));
end
SpectrogramX1Copy = SpectrogramX1;
SpectrogramX1 = 20 * log10(abs(SpectrogramX1));
SpectrogramX1 = SpectrogramX1 - max(SpectrogramX1')';

X = I - 1i * Q;  % Create complex signal

% Ensemble Matrix initialization for complex signal
WindowCollectionMatrixX = zeros(NumberOfTimeWindows, SamplesInWindow * (1 + N));

% Filling the Ensemble Matrix with Zero Padding for complex signal
for kk = 1 : NumberOfTimeWindows   
    TemporaryVector = X(((kk - 1) * SamplesInWindow + 1) : (kk * SamplesInWindow));
    TemporaryZeroVector = zeros(1, SamplesInWindow * N);
    TemporaryVectorZeroPadding = [TemporaryVector, TemporaryZeroVector];
    WindowCollectionMatrixX(kk, :) = TemporaryVectorZeroPadding; 
end

% FFT for complex signal
SpectrogramX2 = zeros(size(WindowCollectionMatrixX));
for kk = 1 : NumberOfTimeWindows
    SpectrogramX2(kk, :) = fft(WindowCollectionMatrixX(kk, :));
end
SpectrogramX2Copy = SpectrogramX2;
% Spectrogram dB conversion for complex signal
SpectrogramX2 = 20 * log10(abs(SpectrogramX2));

% Normalization for complex signal
SpectrogramX2 = SpectrogramX2 - max(SpectrogramX2')';
SpectrogramX2 = CutLowValue(SpectrogramX2,-30);
SpectrogramX1 = CutLowValue(SpectrogramX1,-30);

vdsb = [-fliplr(v) v];  
DirectionSpectrogram = [fliplr(SpectrogramX2) SpectrogramX1];
figure(1)
imagesc(vdsb, t, fliplr(DirectionSpectrogram)); 
axis([-10 10 0 AudioDuration]);
clim([-60 0])
title("Target Velocity CW IF MODE, fc = " + fc / 1e9 + " GHz")
xlabel("Velocity [m/s]", 'FontSize', 12, 'FontWeight', 'bold');
ylabel("Time [s]", 'FontSize', 12, 'FontWeight', 'bold');
hcb = colorbar;
hcb.Title.String = "[dB]";
 
figure(2)
subplot(2,2,1)
plot(t,I); xlabel("time [s]"); ylabel("magnitude"); title("CW LOW IF, REAL PART SIGNAL VS TIME PLOT")
subplot(2,2,2)
plot(t,Q); xlabel("time [s]"); ylabel("magnitude"); title("CW LOW IF, IMAG. PART SIGNAL VS TIME PLOT")
subplot(2,2,3)
plot(t,abs(X)); xlabel("time [s]"); ylabel("magnitude"); title("CW LOW IF, ABS(I+jQ) SIGNAL VS TIME PLOT")
subplot(2,2,4)
plot(t,angle(X)); xlabel("time [s]"); ylabel("magnitude"); title("CW LOW IF, PHASE(I+jQ) SIGNAL VS TIME PLOT")
