function X=pinv_svd(A,tol)
% PINV_SVD	Presudoinverse.
%   X = PINV_SVD(A,tol) Produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD and any
%   singular value less than tolerance 'tol' are treated as zero.
%
%   See also SVD, IMA_SVD, NULL_SVD.

% F. Celi and F. Pasqualetti 2022

[U S V] = svd(A);
r = min(size(S));

while r>0 && abs(S(r,r))<tol
    S(r,r) = 0;
    r = r - 1;
end

Si = inv(S(1:r,1:r));
S(1:r,1:r) = Si;
X = V*S'*U';