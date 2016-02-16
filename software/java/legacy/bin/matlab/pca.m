
function [u1,s1,v1,f] = pca (x, k)
    nrow = size(x,1);
    ncol = size(x,2);

    x = x- repmat(mean(x), nrow,1);
%    x = x./repmat(sqrt(var(x)  + 1/sqrt(nrow)),nrow,1);

    [u,s,v] = svd(x);
    
    f = diag (s);
    f = (f(1:k)'*f(1:k))/(f'*f);
    
    s1 = s (:,1:k);
    v1 = v (:,1:k);
    u1 = u * s1;
    
return