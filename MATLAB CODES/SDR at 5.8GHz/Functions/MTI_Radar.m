function MatOut = MTI_Radar(M)
% Performs Moving Target Indicator of M.
MatrixCopy1 = M( 1:(length(M(:,1))-1) , :);
MatrixCopy2 = M( 2:(length(M(:,1))) , :);

ZeroVec = zeros( length(M(1,:)) , 1 );

MTImatrix = MatrixCopy2-MatrixCopy1;

MatOut =[MTImatrix' ZeroVec]';
end