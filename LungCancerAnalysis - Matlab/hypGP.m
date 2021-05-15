function pvalue = hypGP(N,K,n,k,x)
% N = total
% K = type 1 total
% n = # selected
% k = # type 1 selected
% x = sum to K or n - put numerical value of K or n
warning('off','MATLAB:nchoosek:LargeCoefficient')
pvalue = 0;

for i = k:x
   
    probp = (nchoosek(K,i)*nchoosek((N-K),(n-i)))/(nchoosek(N,n));
    pvalue = pvalue+probp;
    
end
warning('on','MATLAB:nchoosek:LargeCoefficient')

