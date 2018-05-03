function y = mink(A, k)
% returns the smallest k values from a
[~,idx] = sort(A); 
k=min(k,size(idx));
y = idx(1:k);
end