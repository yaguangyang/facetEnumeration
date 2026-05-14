function fOne=firstMinInRow(A)
% All elements of input vector A are { 0, 1, 2 }
% Find the index ii of first A(ii)=1
fMin=100;
rLength=length(A);
for ii=1:rLength
    if A(ii)>=1 & A(ii)<fMin;
        fOne=ii; fMin=A(ii);
    end
end