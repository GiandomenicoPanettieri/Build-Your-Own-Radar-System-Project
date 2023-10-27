%% SAR Range Migration Algorithm (RMA), Group 2 BuildYourOwnRadarSystem
clear all;
close all;
clc;

% constants
c  = 299792458;             % speed of light
SamplingFrequency = 44100;
Tp = 0.02;                  % Upchirp time
fstart =2.408e9;            % starting up-chirp frequency
fstop =2.495e9;             % stopping up-chirp frequency
Deltaf = fstop-fstart;

% Observation time
Trp = 3;                 % segment duration         
SamplesInTrp = Trp*SamplingFrequency;

% Audio manipulation
AudioVector = audioread("SAR_CUSTOM_TEST.m4a");
AudioVector = -(AudioVector)';
BackscateredData = AudioVector(1,:);
Sync = AudioVector(2,:);

% Time axis, for plot purpose
dt = 1/SamplingFrequency;
t_audio = (0 : length(AudioVector(1,:))-1)*dt;
t = 0: dt : Trp-dt;                              % the whole Trp
t_upchirp = 0:dt:Tp-dt;                          % 1 upchirp only

%%
% ------------------------------- PROBLEM ------------------------------- % 
%   We need to run through all the samples. According to the working
%   principle of SAR radar, we move in space the radar, then we turn it on
%   and we make it work as a normal FMCW radar. Then we turn off the FMCW, 
%   we move again the radar and we also turn it back on to reiterate the
%   measure.
%   
%   The problem is:
%   1)  To detect when the FMCW Sync signal is turned on and when it is
%       turned off.
%   2)  To collect the right amount of data from both the Backscatered Data
%       and the Sync. The collected data should be then saved into two
%       different marices.
%
% -------------------------- WORKING PRINCIPLE -------------------------- %
%   We run through all the acquired samples. Two flags are used in this
%   context, namely
%       *flag_wait:
%           =1 indicates that we have collected Trp*SampleFrequency samples
%              so we have to wait for the FMCW to be turned off.
%           =0 indicates that we should have the FMCW turned off since the
%              Sync signal should be flat.
%           /* comment */
%           / problems may occour due to false triggering of flag_wait=0 /
%           / this major issue is caused by spurious samples in the sync /
%           / To solve this problem, secure_counter is added. flag_wait  /
%           / can be put back to 0 only if at least Tp*SamplingFrequency /
%           / non-high and non-low values are counted.                   /
%
% this is done through this custom function
[DataMatrix,SyncMatrix,ParseSyncMatrix,NumberOfUpchirpMatrix] =...
    ExtractSyncedDataSegments(BackscateredData,Sync,Trp,Tp);

figure(1)
segment = 12;
plot(t, SyncMatrix(segment,:)); hold on; plot(t, ParseSyncMatrix(segment,:)); hold on; plot(t,DataMatrix(segment,:))
xlabel("time [s]"); ylabel(" normalized magnitude"); title(" 12th segment example");

% We have to integrate through all the upchirps in a specific position(i.e.
% a row) and put the result of integration into IntegratedDataMatrix (thus we are
% parsing Data matrix)
IntegratedDataMatrix = IntegrateSyncedDataSegments(DataMatrix,ParseSyncMatrix,NumberOfUpchirpMatrix, ...
    Tp,Trp,SamplingFrequency);

%% --------------- Hilbert transform with hann windowing --------------- %%
% In this section we evaluate the Cross-Range matrix i.e. the Hilbert
% transform with Hann windowing of the Data matrix which is a matrix that
% has along the x-axis wave number kr and along y-axis wave number kx.


CrossRangeMatrix = SARCustomHilbert(IntegratedDataMatrix,'hann');

% kr axis
dKr = 4*pi*Deltaf/(c*length(CrossRangeMatrix(1,:)));
Kr = ((4*pi*fstart)/c) : dKr : ((4*pi*fstop)/c)-dKr;

% X axis
lambda = c/fstart;
L = lambda;
dL = L/length(DataMatrix(:,1));
Xa = -L/2 : dL : L/2-dL;

figure(2)
imagesc(Kr, Xa, angle(CrossRangeMatrix)); 
xlabel('Kr [rad/m]'); ylabel('Synthetic Aperture Position Xa[m]'); title('Phase before Column track FFT');
clim([-pi pi]); colorbar; title(colorbar,'[rad]')

%% FFT zero padding of CrossRangeMatrix

DummyZero = zeros(1012,length(CrossRangeMatrix));
CrossRangeMatrix = [DummyZero; CrossRangeMatrix; DummyZero];    %performs zero pdding
CrossRangeMatrixFFT = fftshift(fft(CrossRangeMatrix,[],1),1);
CrossRangeMatrixFFT = CutLowValue(CrossRangeMatrixFFT,-28);

% Kx axis
dKx = 4*pi/(lambda*length(CrossRangeMatrix(:,1)));
Kx = (-2*pi/lambda) : dKx : (2*pi/lambda)-dKx;
% New Kr axis. It has to be adapted with respect to the zero padding
dKr = 4*pi*Deltaf/(c*length(CrossRangeMatrixFFT(1,:)));
Kr = ((4*pi*fstart)/c) : dKr : ((4*pi*fstop)/c)-dKr;

figure(3)
imagesc(Kr,Kx,20*log10(abs(CrossRangeMatrixFFT)));
xlabel('Kr [rad/m]'); ylabel('Kx [rad/m]'); title('Magnitude Cross Range FFT');
clim([-28 14]); colorbar; title(colorbar,'[dB]')

figure(4)
imagesc(Kr,Kx,flipud(angle(CrossRangeMatrixFFT)));
xlabel('Kr [rad/m]'); ylabel('Kx [rad/m]'); title('Angle Cross Range FFT');
clim([-pi pi]); colorbar; title(colorbar,'[rad]')

%% Interpolation
% it is done entirely by means of CustomInterpolation function

[intpolCrossRangeMatrixFFT,Kye] = CustomInterpolation(Kr,Kx,CrossRangeMatrixFFT);

figure(5)
imagesc(Kye,Kx,angle(intpolCrossRangeMatrixFFT));
xlabel('Ky [rad/m]'); ylabel('Kx [rad/m]'); title('Phase After Stolt Interpolation');
clim([-pi pi]); colorbar; title(colorbar,'[rad]')

%% Final Image extraction

c_range_1 = -10;
c_range_2 = 10;
d_range_1 = 1;
d_range_2 = 25;

[FinalTruncatedImage, CrossRange, DownRange] = GetImageFromCustomInterpolation( ...
    intpolCrossRangeMatrixFFT, Kye, c_range_1, c_range_2, d_range_1, d_range_2);

figure(6)
imagesc(CrossRange,DownRange, FinalTruncatedImage); colorbar; clim([-40 0]);
axis([2*c_range_1 2*c_range_2 d_range_1 d_range_2]);
xlabel('CrossRange [m]'); ylabel('DownRange absolute value [m]'); title('Final Image');

