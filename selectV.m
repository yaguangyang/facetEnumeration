function D=selectV(D)
% The elements of the input column vector are {0,1,2,3}
% Replace every 1 by 2 and the rest are unchanged
len=length(D);
for ii=1:len
    if D(ii)==1
        D(ii)=2;
    end
end