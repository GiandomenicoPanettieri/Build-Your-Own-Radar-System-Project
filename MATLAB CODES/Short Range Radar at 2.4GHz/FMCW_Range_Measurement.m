clear all;
close all;
clc;

% constants
c  = 299792458;         % speed of light
fstart =2.408e9;        % starting up-chirp frequency
fstop =2.495e9;         % stopping up-chirp frequency
Deltaf = fstop-fstart;  
Tp = 2e-2;
% Data recovering from audio file
RangeTestFile = audioread("Range_Test_File.m4a");
RangeTestFile = -(RangeTestFile)';                  % adjust the data 
                                                    % taking into account
                                                    % also the -180° shift
                                                    % od the audio card

BackscateredData = RangeTestFile(1,:);
Sync = RangeTestFile(2,:);

% Parsing Sync
ParseSync = -ones(1, length(Sync(1,:)));
for i = 1 : length(Sync(1,:))
    if Sync(i) > 0
        ParseSync(i) = 1;
    elseif Sync(i) < 0
        ParseSync(i) = -1;
    end
end


% Counting the number of upchirps
%
% NumberUpchirps will be used to count the rows of our pectrogram.
% The proposed algorithm runs through all the samples of Sync. If the
% Sync is >0 NumberUpchirps++, we stop adding the NumberUpchirps by setting
% high WaitFor0. When a 0 comes, we set low WaitFor0 and we wait the next
% upchirp
NumberUpchirps = 0;
WaitFor0 = 0;
for i = 1:length(Sync)

    if Sync(i) < 0
        WaitFor0 = 0;
    end
    if Sync(i) > 0 && WaitFor0 == 0
        NumberUpchirps = NumberUpchirps + 1;
        WaitFor0 = 1;
    end

end
% Time axis
SamplingFrequency = 44100;
dt = 1/SamplingFrequency;
AudioDuration = length(Sync)*dt;
t= (0:NumberUpchirps-1)*2*Tp;

% Ensemble Matrix 
Tp = 0.02;          %Upchirp time
EnsembleMatrix = zeros(NumberUpchirps, Tp*SamplingFrequency);
SampleCounter = 1;  %This variable is used to count all BackscateredData samples

for k = 1:NumberUpchirps
    
    %We increase SampleCounter untill Sync goes up 
    while SampleCounter <= length(Sync) && Sync(SampleCounter) <= 0.5
        SampleCounter = SampleCounter + 1;      
    end
    
    % Check if SampleCounter went beyond the length of Sync.
    % This is to avoid errors
    if SampleCounter > length(Sync)
        break; 
    end
    
    % Fill EnsembleMatrix with BackscateredData when Sync is high
    for j = 1 : Tp*SamplingFrequency
        if(SampleCounter+j <= length(ParseSync))
                EnsembleMatrix(k, j) = BackscateredData(SampleCounter + j);
        end
    end
    
    % Move SampleCounter to the start of the next low (<= -0.5) region in Sync
    while SampleCounter <= length(Sync) && Sync(SampleCounter) > -0.5
        SampleCounter = SampleCounter + 1;
    end
end

% Zero Padding
%
% In order to better see the IFFT, a Zero Padding is strongly recomended
N = 20; 
EnsembleMatrixZeroPadding = zeros( NumberUpchirps , (N+1)*length( EnsembleMatrix(1,:) ) );
EnsembleMatrixAddingZero = zeros(1 , N*length(EnsembleMatrix(1,:)));
for i = 1 :NumberUpchirps
    EnsembleMatrixZeroPadding(i,:) = [EnsembleMatrix(i,:) EnsembleMatrixAddingZero];
end

% MS Mean Subtraction
%
% We subtract the mean of each row
for i = 1 : length(EnsembleMatrix(1,:))
    EnsembleMatrixZeroPadding(:,i)=EnsembleMatrixZeroPadding(:,i)-mean(EnsembleMatrixZeroPadding(:,i));
end


%EnsembleMatrixZeroPadding = MTI_Radar(EnsembleMatrixZeroPadding);


% ifft
%
% Here we create the Spectrogram and we turn it into dB
DataSpectrogram = zeros(size(EnsembleMatrixZeroPadding));
for i = 1: NumberUpchirps
    DataSpectrogram(i,:) = 20*log10(abs(ifft(EnsembleMatrixZeroPadding(i,:))))-20*log10(max(abs(ifft(EnsembleMatrixZeroPadding(i,:)))));
end

%DataSpectrogram = DataSpectrogram -max(DataSpectrogram,[],"all");

%cleaning spectrogram
threshold =-50;
indices_below_threshold = find(DataSpectrogram < threshold);

% Replace those elements with the desired value (-60)
DataSpectrogram(indices_below_threshold) = -60;

%range axis
Rmax = c*Tp/2;
dR = c/(2*Deltaf*(N+1));
R = 0:dR:dR*length(DataSpectrogram(1,:))-dR;

figure(1)
imagesc(R,t,DataSpectrogram(:,1:length(R))); clim([threshold 0]); axis([0 100 0 AudioDuration]);
xlabel("Range [m]"); ylabel("time [s]"); colorbar;
title("Range Spectrogram MTI, fstart="+fstart*1e-9+" GHz, fstop="+fstop*1e-9+"GHz");