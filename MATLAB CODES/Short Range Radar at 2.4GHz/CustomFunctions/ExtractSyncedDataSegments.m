function [DataMatrix, SyncMatrix, ParseSyncMatrix, NumberOfUpchirpMatrix] = ExtractSyncedDataSegments(BackscatteredData, Sync, Trp, Tp)
    % EXTRACTSYNCEDDATASEGMENTS Extracts synchronized data segments from backscattered data.
    % [DataMatrix, SyncMatrix, ParseSyncMatrix, NumberOfUpchirpMatrix] = ExtractSyncedDataSegments(BackscatteredData, Sync, Trp, Tp) extracts
    % segments of data from BackscatteredData that are synchronized based on the Sync signal.
    %
    % Inputs:
    % - BackscatteredData: An array containing the backscattered data.
    % - Sync: An array containing the synchronization signal.
    % - Trp: The repetition period of the signal.
    % - Tp: The up-chirp time of the signal.
    %
    % Outputs:
    % - DataMatrix: A matrix where each row is a segment of data from BackscatteredData collected when the synchronization signal was above a certain threshold.
    % - SyncMatrix: A matrix where each row is a segment of the Sync signal corresponding to the data in DataMatrix.
    % - ParseSyncMatrix: Parsed version of SyncMatrix, a square wave signal that goes from -1 to 1 in value.
    % - NumberOfUpchirpMatrix: A column vector where the i-th element is the number of upchirps in the i-th row of DataMatrix.
    %
    % This function works by iterating through all the samples in the backscattered data. It waits until
    % it finds a synchronization signal that is not zero, then starts collecting samples when the signal is above a certain threshold.
    %
    % This function could be useful in applications such as SAR radar signal processing or other areas where
    % you need to synchronize and extract segments of data based on a certain condition.

    SamplingFrequency = 44100;

    SyncMatrix = zeros(1, Trp * SamplingFrequency);
    DataMatrix = zeros(1, Trp * SamplingFrequency);
    CopySyncVector = zeros(1, Trp * SamplingFrequency);
    CopyDataVector = zeros(1, Trp * SamplingFrequency);
    flag_wait = 1;
    flag_start = 0;
    secure_counter0 = 0;
    secure_counter1 = 0;

    for i = 1:length(BackscatteredData(1,:))  % Iterate through all the samples

        if (Sync(i) < 0.1) && (Sync(i) > -0.1) && flag_wait == 1 && flag_start == 0
            secure_counter0 = secure_counter0 + 1;
            if secure_counter0 == Trp * SamplingFrequency
                flag_wait = 0;
            end
        elseif (Sync(i) > 0.8) && flag_wait == 0 && flag_start == 0
            secure_counter1 = secure_counter1 + 1;
            if secure_counter1 == Trp * SamplingFrequency
                flag_start = 1;
            end
        end

        if flag_start == 1
            flag_start = 0;
            for jj = 1:Trp * SamplingFrequency
                CopySyncVector(1, jj) = Sync(1, jj + i);
                CopyDataVector(1, jj) = BackscatteredData(1, jj + i);
            end
            SyncMatrix = [SyncMatrix; CopySyncVector];
            DataMatrix = [DataMatrix; CopyDataVector];
            i = i + jj;  % Update i to the correct value
            secure_counter0 = 0;
            flag_wait = 1;
            secure_counter1 = 0;
        end
    end

    % Remove the first row of zeros from both SyncMatrix and DataMatrix
    SyncMatrix = SyncMatrix(2:end, :);
    DataMatrix = DataMatrix(2:end, :);

    % Parse all syncs
    ParseSyncMatrix = -ones(size(SyncMatrix));
    NumberOfUpchirpMatrix = zeros(length(SyncMatrix(:, 1)), 1);
    for jj = 1:length(SyncMatrix(:, 1))
        for i = 1:Trp * SamplingFrequency
            if SyncMatrix(jj, i) > 0.5
                ParseSyncMatrix(jj, i) = 1;
            elseif SyncMatrix(i) < -0.5
                ParseSyncMatrix(jj, i) = -1;
            end
        end
        NumberOfUpchirpMatrix(jj) = UpchirpsCount(SyncMatrix(jj, :));  % Collect the number of upchirps into this matrix
    end

end
