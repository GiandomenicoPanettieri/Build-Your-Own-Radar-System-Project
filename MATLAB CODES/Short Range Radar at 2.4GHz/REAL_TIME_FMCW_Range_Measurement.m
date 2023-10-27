clear all;
close all;
clc;

% constants
c  = 299792458;         % speed of light
fstart =2.408e9;        % starting up-chirp frequency
fstop =2.495e9;         % stopping up-chirp frequency
Deltaf = fstop-fstart;  

% Time axis
SamplingFrequency = 44100;
Tp = 0.02;                                  % Upchirp time
dt = 2*Tp;
N = 4;
SpectrogramTimeWindow= 10;                  % spectrogram seconds
t = (0:SpectrogramTimeWindow/(2*Tp))*dt;


% Spectrogram initialization
TotalSpectrogram = zeros(SpectrogramTimeWindow/(2*Tp), (N+1)*Tp*SamplingFrequency);

% range axis
Rmax = c*Tp/2;
dR = c/(2*Deltaf*(N+1));
RefreshTime=0.36;                                    % c RefreshTime of the Radar
TotalRows = SpectrogramTimeWindow/(2*Tp);           % total rows of the ensemble matrix
R = 0:dR:dR*length(TotalSpectrogram(1,:))-dR;

%THIS SIMULATES A REAL-TIME MESUREMENT. YOU CAN READ Tp*SamplingFrequency
%samples from USB port as well!
k=1;
while 1 
    RangeTest = audioread('Range_Test_File.m4a', [(RefreshTime*SamplingFrequency*k +1) ...
        (RefreshTime*SamplingFrequency*(k+1))]);
    RangeTest = -(RangeTest)';                      % adjust the data 
                                                    % taking into account

                                                    % od the audio card

    BackscateredData = RangeTest(1,:);
    Sync = RangeTest(2,:);

     % Parsing Sync
    ParseSync = -ones(1, length(Sync(1,:)));
    for i = 1 : length(Sync(1,:))
        if Sync(i) > 0
            ParseSync(i) = 1;
        elseif Sync(i) < 0
            ParseSync(i) = -1;
        end
    end
    
    NumberUpchirps = UpchirpsCount(ParseSync);

    LocalEnsambleMatrix = EnsembleMatrixZPfill(ParseSync,BackscateredData,NumberUpchirps,N);
    for i = 1 : length(LocalEnsambleMatrix(1,:))
        LocalEnsambleMatrix(:,i)=LocalEnsambleMatrix(:,i)-mean(LocalEnsambleMatrix(:,i));
    end
    %LocalEnsambleMatrix = MTI_Radar(LocalEnsambleMatrix);
    LocalSpectrogram = zeros(size(LocalEnsambleMatrix));
    for i = 1 : NumberUpchirps
        LocalSpectrogram(i,:) = abs( ifft( LocalEnsambleMatrix(i,:) ) )/max( abs( ifft( LocalEnsambleMatrix(i,:) ) ) );
    end
    LocalSpectrogram = 20*log10(LocalSpectrogram);
    LocalSpectrogram = LocalSpectrogram -max(LocalSpectrogram,[],"all");

    TotalSpectrogram = PutOnTop(TotalSpectrogram,LocalSpectrogram);

    TotalSpectrogram = CutLowValue(TotalSpectrogram,-20);
    
    imagesc(R,t,TotalSpectrogram(:,1:length(R))); clim([-50 0]); axis([0 100 0 t(length(t))]);
    xlabel("Range [m]"); ylabel("time [s]");
    title("Range Spectrogram, fstart="+fstart*1e-9+" GHz, fstop="+fstop*1e-9+"GHz");
    pause(RefreshTime/5)

    k=k+1;
end
