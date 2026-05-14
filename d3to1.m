function D=d3to1(D)
% All elements of input matrix D are either 0, 1, 2, 3
% This function changes all elements of 3 to 1, the rest elements are
% unchanged
[m,n]=size(D);
for ii=1:m
    for jj=1:n
        if D(ii,jj)==3
            D(ii,jj)=1;
        end
    end
end