function IntegratedDataMatrix = IntegrateSyncedDataSegments(DataMatrix,ParseSyncMatrix,NumberOfUpchirpMatrix,Tp,Trp,SamplingFrequency)
%INTEGRATESYNCEDDATASEGMENTS Integrates, through all the upchirps, each row of a matrix DataMatrix that have synced data segments 
%   organized through each row.
%   This function takes in a DataMatrix, a ParseSyncMatrix, a matrix containing the number of upchirps,
%   and parameters Tp, Trp, and SamplingFrequency. It then integrates through all the upchirps in a specific
%   position (i.e., a row) and puts the result of integration into IntegratedDataMatrix (thus parsing Data matrix).
%
%   Inputs:
%       DataMatrix - A matrix which have segments of BackscateredData organized on each
%                    row. See ExtractSyncedDataSegments.
%       ParseSyncMatrix - The parse sync matrix. It must be the pase
%       version!
%       NumberOfUpchirpMatrix - The matrix containing the number of upchirps.
%       Tp - Parameter related to the duration of the upchirp.
%       Trp - Parameter related to the segment duration.
%       SamplingFrequency - The sampling frequency.
%
%   Outputs:
%       IntegratedDataMatrix - The updated data matrix after integration through all the upchirps.
%


% We have to integrate through all the upchirps in a specific position(i.e.
% a column) and put the result of integration into DataMatrix (thus we are
% parsing Data matrix)
    DataMatrixCopy = DataMatrix;                        % we copy Data matrix
    DataMatrix = zeros(length(NumberOfUpchirpMatrix(:,1)),Tp*SamplingFrequency);        % new expected 
                                                        % dimension of
                                                        % DataMatrix
    UpchirpDataMatrix = zeros(1,Trp*SamplingFrequency);
    SingleUpchirpDataVector = zeros(1,Tp*SamplingFrequency);
    AddingUpchirpDataVector = zeros(1,Tp*SamplingFrequency);
    
    for jj = 1 : length(NumberOfUpchirpMatrix(:,1)) %running through all the upchirps
        % Custom function EnsembleMatrixZPfill provides a very easy way to
        % solve the problem. The ensemble matrix is the matrix that contains in
        % its rows the samples of the backscatered data when the sync is high
        LocalEnsembleMatrix = EnsembleMatrixZPfill(ParseSyncMatrix(jj,:), ...
            DataMatrixCopy(jj,:),NumberOfUpchirpMatrix(jj,:),0);
        for kk = 1 : NumberOfUpchirpMatrix(jj,:)-1
            DataMatrix(jj,:) =DataMatrix(jj,:)+LocalEnsembleMatrix(kk,:);
        end
        DataMatrix(jj,:)=DataMatrix(jj,:)/(NumberOfUpchirpMatrix(jj)-1);
    end
    IntegratedDataMatrix = DataMatrix;
end