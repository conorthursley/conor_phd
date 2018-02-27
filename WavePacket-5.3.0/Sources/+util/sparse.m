function [ Y ] = sparse( X )
% Sparsify an approximately sparse matrix by squeezing out any 
% elements whose absolute values are below (relative) threshold. 
% If X is complex then this criterion is applied to real and  
% imaginary parts separately.
% If the resulting density is below 10% the returned matrix will
% be in sparse form; otherwise a full matrix will be returned.

tol = eps('single');

if isreal(all(X(:)))
    thresh = tol * max(abs(X(:)));
    Y = X.*(abs(X)>thresh);
else
    thresh1 = tol * max(abs(real(X(:))));
    thresh2 = tol * max(abs(imag(X(:))));
    Y = ...
        real(X).*(abs(real(X))>thresh1) + ...
        imag(X).*(abs(imag(X))>thresh2) * 1i;
end

if nnz(Y)/numel(Y)<0.1
    Y = sparse(Y);
end

% util.disp (['Density of matrix before: ' num2str(nnz(X)/numel(X))])
% util.disp (['Density of matrix after:  ' num2str(nnz(Y)/numel(Y))])




