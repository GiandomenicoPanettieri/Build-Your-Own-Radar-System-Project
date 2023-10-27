function Mout= PutOnTop(M,V)
%   Mout= PutOnTop(M,V) Puts matrix V on top of matrix M without changing the size of M
 
    M = [V ; M ];
    M = M(1:end-length(V(:,1)),:);
Mout=M;
end