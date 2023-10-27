function [intpolCrossRangeMatrixFFT,Kye] = CustomInterpolation(Kr,Kx,CrossRangeMatrixFFT)
%CUSTOMINTERPOLATION computes the interpolated Cross-Range Matrix FFT and returns the 
% corresponding cross-range wave numbers.
%
% Syntax:
%   [intpolCrossRangeMatrixFFT, Kye] = KyMatrix(Kr, Kx, CrossRangeMatrixFFT)
%
% Input Arguments:
%   - Kr: Vector of range wave numbers.
%   - Kx: Vector of cross-range wave numbers.
%   - CrossRangeMatrixFFT: Original Cross-Range Matrix in the wave number domain.
%
% Output Arguments:
%   - intpolCrossRangeMatrixFFT: Interpolated Cross-Range Matrix in the wave number domain.
%   - Kye: Vector of interpolated cross-range wave numbers.
%
% Description:
%   This function takes the range wave numbers (Kr), cross-range wave numbers (Kx), and
%   the original Cross-Range Matrix in the wave number domain (CrossRangeMatrixFFT) as
%   input. It computes the cross-range wave numbers (Kyr) for each Kx value and then
%   interpolates the Cross-Range Matrix to obtain intpolCrossRangeMatrixFFT using the
%   interpolated Kyr values.
%
%   NaN values in the interpolated Cross-Range Matrix are replaced with a small
%   positive value (1e-30) to avoid numerical issues.
%
%   The function returns the interpolated Cross-Range Matrix (intpolCrossRangeMatrixFFT)
%   and the corresponding interpolated cross-range wave numbers (Kye).

    Ky = zeros( length(Kx), length(Kr)  );
    for jj = 1 : length(Kx)
        Ky(jj,:) = sqrt(Kr.^2-Kx(jj).^2);
    end

    Kye_min =min(Ky,[],'all');
    Kye_max =max(Ky,[],'all');
    dKye = (Kye_max-Kye_min)/(length(Kx)/2);
    Kye = Kye_min : dKye : Kye_max-dKye;
    
    % interpolated CrossRangeMatrixFFT
    intpolCrossRangeMatrixFFT = zeros(length(CrossRangeMatrixFFT(:,1)),length(Kye));
    for jj = 1 : length(CrossRangeMatrixFFT(:,1))
        intpolCrossRangeMatrixFFT(jj,:) = interp1(Ky(jj,:),CrossRangeMatrixFFT(jj,:),Kye);
    end
    WhereNaN = isnan(intpolCrossRangeMatrixFFT);
    intpolCrossRangeMatrixFFT(WhereNaN)=1e-30;
end