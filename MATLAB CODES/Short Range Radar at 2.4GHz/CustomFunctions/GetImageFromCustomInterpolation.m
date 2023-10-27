function [FinalTruncatedImage, CrossRange, DownRange] = GetImageFromCustomInterpolation( ...
    intpolCrossRangeMatrixFFT, Kye, c_range_1, c_range_2, d_range_1, d_range_2)
%GETIMAGEFROMCUSTOMINTERPOLATION Process radar data and generate a truncated image
%   This function takes radar data in the form of a 2D FFT (Fast Fourier Transform)
%   matrix, performs various transformations and truncations, and generates a
%   final truncated radar image.
%
%   Input Arguments:
%   - intpolCrossRangeMatrixFFT: 2D FFT matrix representing the radar data.
%   - Kye: Cross-range wavenumber vector.
%   - c_range_1, c_range_2: Cross-range limits for truncation.
%   - d_range_1, d_range_2: Down-range limits for truncation.
%
%   Output Arguments:
%   - FinalTruncatedImage: The final truncated radar image.
%   - CrossRange: Vector representing cross-range positions in the image.
%   - DownRange: Vector representing down-range positions in the image.
%
%   Detailed Explanation:
%   - The function starts by performing a 2D IFFT (Inverse Fast Fourier Transform)
%     on the input radar data matrix, with zero-padding.
%   - Flips the resulting matrix to adjust for orientation.
%   - The truncated image region is defined based on specified cross-range and
%     down-range limits.
%   - Cross-range and down-range vectors are generated based on the specified limits.
%   - The final truncated radar image is computed by applying specific operations
%     on the truncated data, including multiplication by the squared absolute
%     values of the down-range positions and logarithmic scaling.


    % constants
    c  = 299792458;             % speed of light
    SamplingFrequency = 44100;
    Tp = 0.02;                  % Upchirp time
    fstart =2.408e9;            % starting up-chirp frequency
    fstop =2.495e9;             % stopping up-chirp frequency
    Deltaf = fstop-fstart;
    lambda = c/fstart;

    % Range Data Matrix will be our final matrix
    % Perform ifft2D. Zero Padding is set to 4 by default.
    RangeDataMatrix = ifft2(intpolCrossRangeMatrixFFT, ...
        4*length(intpolCrossRangeMatrixFFT(:,1)), ...
        4*length(intpolCrossRangeMatrixFFT(1,:)) );
    
    % Rotate -90° RangeDataMatrix and flip left right
    RangeDataMatrix = fliplr(RangeDataMatrix);
    RangeDataMatrix = rot90(rot90(rot90(RangeDataMatrix, 1)));
    
    % definition of the truncated Image
    dx = lambda/2;                        % displacement of the radar on top of the rail (OK)
    dfy = (c/(2*pi))*(Kye(end)-Kye(1)+18);     
    Rmax = (c/(2*dfy))*(length(RangeDataMatrix(:,1)));
    Rail_Rmax = length(intpolCrossRangeMatrixFFT(:,1))*dx;

    d_index_1 = round((size(RangeDataMatrix,1)/Rmax)*d_range_1);
    d_index_2 = round((size(RangeDataMatrix,1)/Rmax)*d_range_2);
    c_index_1 = round((size(RangeDataMatrix,2)/Rail_Rmax)*(c_range_1+(Rail_Rmax/2)));
    c_index_2 = round((size(RangeDataMatrix,2)/Rail_Rmax)*(c_range_2+(Rail_Rmax/2)));
    
    TruncatedRangeDataMatrix = RangeDataMatrix(d_index_1:d_index_2,c_index_1:c_index_2);
    TruncatedRangeDataMatrix = fliplr(TruncatedRangeDataMatrix);
    
    % cross range vector 
    dcr = (c_range_2 - c_range_1)/(c_index_2-c_index_1);
    CrossRange = c_range_1 : dcr : c_range_2;
    
    % down range vector 
    ddr = (d_range_2 - d_range_1)/(d_index_2-d_index_1);
    DownRange = d_range_1 : ddr : d_range_2;
    
    FinalTruncatedImage = zeros(length(DownRange),length(CrossRange));
    for ii = 1 : length(TruncatedRangeDataMatrix(1,:))
        FinalTruncatedImage(:,ii) = TruncatedRangeDataMatrix(:,ii) .* (abs(DownRange')).^2;
    end
    FinalTruncatedImage=20*log10(abs(FinalTruncatedImage)/max(abs(FinalTruncatedImage),[],'all'));
end