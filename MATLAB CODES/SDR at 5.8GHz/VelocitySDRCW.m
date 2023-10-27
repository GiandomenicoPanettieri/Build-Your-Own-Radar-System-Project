clear all;
close all;
clc;

%Constants. Suggested parameters 
fc = 5.8e9;
Tp = 0.1;
N = 4;
c  = 299792458;                     % speed of light

SamplingFrequency=44100;            % Referred to the audio signal

%%                                           REAL PART
% Data from Audio
I = audioread("RealPartCW.wav");
I = I';

% mean subtraction (MS)
I = I-mean(I);
AudioDuration = length(I)/SamplingFrequency;  %time duration of input data
% time axis
dt = 1/SamplingFrequency;
t = 0 : dt : AudioDuration-dt;
%frequency and velocity axis
BandWidth = SamplingFrequency;
bin = SamplingFrequency/(length(I)*(N+1)); 
f = 0 : bin : BandWidth-bin;                        
v = f * (c/(2*fc));                                 
        
% Ensemble Matrix initialisation
NumberOfTimeWindows = floor(AudioDuration/(Tp));

SamplesInWindow = Tp*SamplingFrequency;         
WindowCollectionMatrixI = zeros(NumberOfTimeWindows, SamplesInWindow*(1+N));

% Filling the Ensemble Matrix With also Zero Padding

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
% SpectrogramI dB conversion
 SpectrogramI = 20*log10( abs(SpectrogramI) );
 SpectrogramI = CutLowValue(SpectrogramI,-60);
% Normalization 2 as default

SpectrogramI = SpectrogramI-max(SpectrogramI')';
figure(1)
imagesc(v,t,SpectrogramI); axis([0 5 0 AudioDuration]); clim([-10 0])
text = "SDR Real Part Single Target CW Low IF, WindowTime Tp = " +Tp+", CarrierFrequency fc="+fc/1e9+" GHz";
title(text)        
        
% This is just plot fashion.
xlabel("Velocity [m/s]",'FontSize',12,'FontWeight','bold');
ylabel("Time [s]",'FontSize',12,'FontWeight','bold');
hcb=colorbar;
hcb.Title.String = "[dB]";

%%                                           IMAG PART
% Data from Audio
Q = audioread("ImaginaryPartCW.wav");
Q = Q';

% mean subtraction (MS)
Q = Q-mean(Q);
AudioDuration = length(Q)/SamplingFrequency;  %time duration of input data
% time axis
dt = 1/SamplingFrequency;
t = 0 : dt : AudioDuration-dt;
%frequency and velocity axis
BandWidth = SamplingFrequency;
bin = SamplingFrequency/(length(Q)*(N+1)); 
f = 0 : bin : BandWidth-bin;                        
v = f * (c/(2*fc));                                 
        
% Ensemble Matrix initialisation
NumberOfTimeWindows = floor(AudioDuration/(Tp));

SamplesInWindow = Tp*SamplingFrequency;         
WindowCollectionMatrixQ = zeros(NumberOfTimeWindows, SamplesInWindow*(1+N));

% Filling the Ensemble Matrix With also Zero Padding

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
% SpectrogramI dB conversion
 SpectrogramQ = 20*log10( abs(SpectrogramQ) );
% Normalization 2 as default

SpectrogramQ = SpectrogramQ-max(SpectrogramQ')';
SpectrogramQ = CutLowValue(SpectrogramQ,-10);
figure(2)
imagesc(v,t,SpectrogramQ); axis([0 5 0 AudioDuration]); clim([-10 0])
text = "SDR Imaginary Part Single Target CW, WindowTime Tp = " +Tp+", CarrierFrequency fc="+fc/1e9+" GHz";
title(text)        

% This is just plot fashion.
xlabel("Velocity [m/s]",'FontSize',12,'FontWeight','bold');
ylabel("Time [s]",'FontSize',12,'FontWeight','bold');
hcb=colorbar;
hcb.Title.String = "[dB]";

%% Spectrogram of Complex Signal
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
% Spectrogram dB conversion for complex signal
SpectrogramX1 = 20 * log10(abs(SpectrogramX1));

% Normalization for complex signal
SpectrogramX1 = SpectrogramX1 - max(SpectrogramX1')';
%% Spectrogram of Complex Signal
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
SpectrogramX2 = CutLowValue(SpectrogramX2,-20);
SpectrogramX1 = CutLowValue(SpectrogramX1,-20);

vdsb = [-fliplr(v) v];
DirectionSpectrogram = [fliplr(SpectrogramX2) SpectrogramX1];
figure(5)
imagesc(vdsb, t, DirectionSpectrogram); 
axis([-5 5 0 AudioDuration]);
clim([-30 0])
title("SDR Velocitiy measurement in CW mode")
xlabel("Velocity [m/s]", 'FontSize', 12, 'FontWeight', 'bold');
ylabel("Time [s]", 'FontSize', 12, 'FontWeight', 'bold');
hcb = colorbar;
hcb.Title.String = "[dB]";
 
figure(6)
subplot(2,2,1)
plot(t,I); xlabel("time [s]"); ylabel("magnitude"); title("CW, REAL PART SIGNAL VELOCITY VS TIME PLOT")
subplot(2,2,2)
plot(t,Q); xlabel("time [s]"); ylabel("magnitude"); title("CW, IMAG. PART SIGNAL VELOCITY VS TIME PLOT")
subplot(2,2,3)
plot(t,abs(X)); xlabel("time [s]"); ylabel("magnitude"); title("CW, ABS(I+jQ) SIGNAL VELOCITY VS TIME PLOT")
subplot(2,2,4)
plot(t,angle(X)); xlabel("time [s]"); ylabel("magnitude"); title("CW, PHASE(I+jQ) SIGNAL VELOCITY VS TIME PLOT")
