function NumberUpchirpsOut = UpchirpsCount(Sync)
% Calculates the Number of UpChirps of Sync.
% Sync is expected to be a row vector. The input vector mus swing from -1
% to 1.
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
    
    % Counts the number of upchirps
    NumberUpchirpsOut = NumberUpchirps;

end