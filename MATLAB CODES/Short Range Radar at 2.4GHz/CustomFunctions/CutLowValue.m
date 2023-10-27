function B = CutLowValue(A,threshold)
%CUTLOWVALUE function replaces values in A that are less than or equal to a specified threshold with a very small value (1e-30).
% Inputs:
%   A - Input matrix.
%   threshold - Threshold for cutting low values.
% Output:
%   B - Output matrix with low values replaced.
%
    find_in_A = find(A<=threshold); 
    A(find_in_A) = -120;
    B = A;
end