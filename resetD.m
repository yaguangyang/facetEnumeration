function D=resetD(D)
[m,n]=size(D);
for ii=1:m
    for jj=1:n
        if D(ii,jj)==2
            D(ii,jj)=1;
        end
    end
end