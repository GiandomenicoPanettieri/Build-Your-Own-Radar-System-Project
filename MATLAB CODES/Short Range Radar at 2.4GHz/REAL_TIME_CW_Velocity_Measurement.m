clear all;
close all;
clc;

% Constants
fc = 2.43e9;
Tp = 0.1;                           % Window
N = 4;                              % Zero Padding
c  = 299792458;                     % speed of light
SamplingFrequency=44100;            % Referred to the audio signal

% time axis
dt = 1/SamplingFrequency;
SpectrogramTimeWindow= 5;
t = 0 : dt : SpectrogramTimeWindow-dt; 
SamplesInWindow = Tp*SamplingFrequency;

%frequency and velocity
BandWidth = SamplingFrequency;
bin = SamplingFrequency/(SamplesInWindow*(1+N));    % bin is the Sampling Frequency
                                                    % over the total number
                                                    % of samples after zero
                                                    % padding!

f = 0 : bin : BandWidth-bin;                        % Frequency axis
v = f * (c/(2*fc));                                 % Velocity axis

        

NumberOfTimeWindows = floor(SpectrogramTimeWindow/(Tp));    % Number of time windows that
                                                            % we can construct from our
                                                            % signal
                                                          
Spectrogram_out = (-50)*ones(NumberOfTimeWindows, SamplesInWindow*(1+N));
TemporaryZeroVector = zeros(1, SamplesInWindow*N);

k=0;

%THIS SIMULATES A REAL-TIME MESUREMENT. YOU CAN READ Tp*SamplingFrequency
%samples from USB port as well!
NameOfAudio = 'Velocity_Test_File.m4a';

while 1
  % Signal extraction and zero padding
    VelocityTest = audioread(NameOfAudio, [(Tp*SamplingFrequency*k +1) ...
        (Tp*SamplingFrequency*(k+1))]);
    VelocityTest = -(VelocityTest)';
    
    
    VelocityData = VelocityTest(1,:);
    VelocityDataPadded = [ VelocityData TemporaryZeroVector ];
    
  % performs the FFT data by data and directly normalizes it
    NewSpectrogrmRow = 20*log10(abs(fft(VelocityDataPadded))/max(abs(fft(VelocityDataPadded))));

  % puts on top of spectrogram the new data
    Spectrogram_out = PutOnTop(Spectrogram_out,NewSpectrogrmRow); %help PutOnTop
    Spectrogram_out = CutLowValue(Spectrogram_out,-20);

    figure(1)
    imagesc(v,t,Spectrogram_out); clim([-50 -10]);
    axis([ v(1) 30 0 SpectrogramTimeWindow])
    xlabel("Velocity [m/s]"); ylabel("time [s]");
    title("Velocity Spectrogram");

    figure(2)
    plot(v,NewSpectrogrmRow);axis([ 0 30 -50 0]);
    xlabel("Velocity [m/s]"); ylabel("dB");
    title("Velocity Spectrogram");
    pause(Tp)
    k= k + 1;
end

