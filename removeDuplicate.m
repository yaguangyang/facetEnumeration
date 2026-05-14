function A=removeDuplicate(A)
Ac = zeros(size(A));
for k1 = 1:length(A)-1
    Ac(k1) = sum(A(k1) == A(k1+1:end));
end
out = A(Ac == 0);
A = out;