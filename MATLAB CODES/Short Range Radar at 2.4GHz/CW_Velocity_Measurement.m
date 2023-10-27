clear all;
close all;
clc;

%Constants. Suggested parameters 
fc = 5.8e9;
Tp = 0.1;
N = 4;
c  = 299792458;                     % speed of light

SamplingFrequency=44100;            % Referred to the audio signal

% Data from Audio
% 
% An audio has Stereo/Mono options. In this case we have a stereo file,
% thus the audioread() function provides 2 column vectors. Only one of these two
% vectors contains our data. It turns out that it is column 1
AudioVector = audioread("Velocity_Test_File.m4a");
AudioVector = -AudioVector';    % In this way we take care of the 
                                % 180° phase-shifted signal coming from
                                % the soundcard and we transpose the 
                                % AudioVector in order to have datas 
                                % displayed as row vectors

AudioVector = AudioVector(1,:); % We only select datas from the first row
                                % since the second is only noise

% mean subtraction (MS)
AudioVector = AudioVector-mean(AudioVector);    % Mean subtraction for 
                                                % the entire AudioVector
                                
AudioDuration = length(AudioVector)/SamplingFrequency;  %time duration of input data


% time axis
dt = 1/SamplingFrequency;
t = 0 : dt : AudioDuration-dt;

% plot(t,AudioVector)           %uncomment to see the signal in time 


%frequency and velocity axis
BandWidth = SamplingFrequency;
bin = SamplingFrequency/(length(AudioVector)*(N+1)); % bin is the Sampling Frequency
                                                    % over the total number
                                                    % of samples after zero
                                                    % padding!
f = 0 : bin : BandWidth-bin;                        % Frequency axis
v = f * (c/(2*fc));                                 % Velocity axis
        

%%
% Ensemble Matrix initialisation
NumberOfTimeWindows = floor(AudioDuration/(Tp));% Number of time windows that
                                                % we can construct from our
                                                % signal

SamplesInWindow = Tp*SamplingFrequency;         % How many samples we have in 
                                                % a window. This should
                                                % always be an integer, so
                                                % if Tp<1 DO NOT PROVIDE TP WITH
                                                % MORE THAN 3 DIGITS!
WindowCollectionMatrix = zeros(NumberOfTimeWindows, SamplesInWindow*(1+N));

% Filling the Ensemble Matrix With also Zero Padding
%
% This is the way we perform Zero Padding in this script
% If A = [1 2 3] and we want to add a certain chunck of zeros
% ZeroVector = zeros(1, N*length(A));
% AZeroPadding = [A ZeroVector];
%
% We have to do this procedure for all rows of our EnsembleMatrix

for kk = 1 : NumberOfTimeWindows
    
    % TemporaryVector is a local utilised vector, it contains
    % the AudioVector datas associated to the kk-th Window.
    TemporaryVector = AudioVector( ((kk-1)*SamplesInWindow +1) : ((kk)*SamplesInWindow) );
    %TemporaryVector = TemporaryVector - mean(TemporaryVector);  % MeanSubtraction Window per Window
    TemporaryZeroVector = zeros(1, N*length(TemporaryVector));
    TemporaryVectorZeroPadding = [TemporaryVector TemporaryZeroVector];
    
    WindowCollectionMatrix(kk,:) = TemporaryVectorZeroPadding; 
end

%%
% FFT 
Spectrogram = zeros(size(WindowCollectionMatrix));  % According to the FFT 
                                                    % algorithm, Spectrogram   
                                                    % has the same size of
                                                    % WindowCollectionMatrix

for kk = 1 : NumberOfTimeWindows
    
    Spectrogram(kk,:) = fft(WindowCollectionMatrix(kk,:));

end

% Spectrogram dB conversion

 Spectrogram = 20*log10( abs(Spectrogram) );



% Normalization 2 as default
%
% This time we normalise row by row dividing all the row by the 
% maximum number in that row
Spectrogram = Spectrogram-max(Spectrogram')';
figure(2)
imagesc(v,t,Spectrogram); axis([0 14 2 AudioDuration]); clim([-10 0])
text = "CW Velocity Measurement, WindowDuration" +Tp+", CenterFrequency fc="+fc/1e9+" GHz";
title(text)        
        

% This is just plot fashion.
xlabel("Velocity [m/s]",'FontSize',12,'FontWeight','bold');
ylabel("Time [s]",'FontSize',12,'FontWeight','bold');
hcb=colorbar;
hcb.Title.String = "[dB]";
