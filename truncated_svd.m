function [U, S, V] = truncated_svd(A, eps)
    %obtain all the singular values of A
    sigma = svd(A);
    
    %obtain the number of the singular values>=eps of A
    sigma_large = sigma(sigma>=eps);
    len = length(sigma_large);
    
    %calculate the trancated svd of A
    [U,S,V] = svds(A,len);
end