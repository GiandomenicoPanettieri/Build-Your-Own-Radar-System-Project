function outmatrix = EnsembleMatrixZPfill(ParseSync,BackscateredData,NumberUpchirps,N)
% Creates the EnsembleMatrix from BackscateredData, with zero padding specified by N and using
% ParseSync.
% N=0 means No zero padding
% ParseSync should go from -1 to 1. Error occurs if it goes from 0 to 1
    Tp = 0.02;
    SamplingFrequency = 44100;
    LocalEnsembleMatrix = zeros(NumberUpchirps, Tp*SamplingFrequency);
    SampleCounter = 1;  %This variable is used to count all BackscateredData samples
    
    for kk = 1:NumberUpchirps
        
        %We increase SampleCounter untill Sync goes up 
        while SampleCounter <= length(ParseSync) && ParseSync(SampleCounter) <= 0.5
            SampleCounter = SampleCounter + 1;      
        end
        
        % Check if SampleCounter went beyond the length of Sync.
        % This is to avoid errors
        if SampleCounter >= length(ParseSync)
            break; 
        end
        
        % Fill EnsembleMatrix with BackscateredData when Sync is high
        for j =1 : Tp*SamplingFrequency
            if(SampleCounter+j <= length(ParseSync))
                LocalEnsembleMatrix(kk, j) = BackscateredData(SampleCounter + j);
            end
        end
        
        % Move SampleCounter to the start of the next low (<= -0.5) region in Sync
        while SampleCounter < length(ParseSync) && ParseSync(SampleCounter) > -0.5
            SampleCounter = SampleCounter + 1;
        end
    end
    
    % Zero Padding
    %
    % In order to better see the IFFT, a Zero Padding is strongly recomended
    % N = 4; 
    LocalEnsembleMatrixZeroPadding = zeros( NumberUpchirps , (N+1)*length( LocalEnsembleMatrix(1,:) ) );
    LocalEnsembleMatrixAddingZero = zeros(1 , N*length( LocalEnsembleMatrix(1,:) ) );
    for i = 1 :NumberUpchirps
        LocalEnsembleMatrixZeroPadding(i,:) = [LocalEnsembleMatrix(i,:) LocalEnsembleMatrixAddingZero];
    end
    
    outmatrix = LocalEnsembleMatrixZeroPadding ;
    
end