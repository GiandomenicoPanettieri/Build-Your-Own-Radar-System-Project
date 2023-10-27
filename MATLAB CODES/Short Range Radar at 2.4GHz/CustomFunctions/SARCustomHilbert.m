function A_Hilbert_out = SARCustomHilbert(A,hann_option)
%SARCustomHilbert performs the Hilbert transform on the input matrix A.
% The function can also applie a Hann window if the 'hann_option' is set to 'hann'.
%
% Inputs:
%   A - Input matrix. 
%   hann_option - Optional argument. If set to 'hann', a Hann window is applied to the result of the Hilbert transform
%   row by row          
% Outputs:
%   A_Hilbert_out - The output of the Hilbert transform. If 'hann_option' is set to 'hann', this will be the windowed result.
%

    
    if nargin < 1  
        error ('A is a required input') 
    end  
    if nargin < 2  
      hann_option = 'a';
    end

    % Hilbert transform
    % 1 FFT on each row (position)
    A_FFT = zeros(size(A));
    for jj = 1:length(A(:,1))
        A_FFT(jj,:) = fft(A(jj,:));
    end
    A_FFT_positive_frequency = ...
        A_FFT(:,1:length(A(1,:))/2);
    % 2 IFFT on positive frequency of each row
    A_Hilbert = zeros(size(A_FFT_positive_frequency));
    for jj = 1:length(A(:,1))
        A_Hilbert(jj,:) = ifft(A_FFT_positive_frequency(jj,:));
    end
    
    % Define the value to replace NaN with
    replacementValue = 1e-30;
    
    % Create a logical mask for NaN values in A
    nanMask = isnan(A_Hilbert);
    
    % Replace NaN values with the replacementValue
    A_Hilbert(nanMask) = replacementValue;

    if hann_option == 'hann'
        % Hann windowing row by row 
        window = hann(size(A_Hilbert,2));
        A_windowed_Hilbert = zeros(size(A_Hilbert));
        for i = 1:size(A, 1)
            A_windowed_Hilbert(i, :) = A_Hilbert(i, :) .* window';
        end
        A_Hilbert_out = A_windowed_Hilbert;
    else
        A_Hilbert_out = A_Hilbert;
    end
end