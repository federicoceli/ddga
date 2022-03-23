function Z = null_svd(A,tol)
% NULL_SVD	Null space.
%   Z = NULL_SVD(A,tol) is an orthonormal basis for the null space of A.
%   That is,  A*Z has negligible elements, size(Z,2) is the nullity of A, 
%   and Z'*Z = I. The computation is based on SVD and any
%   singular value less than tolerance 'tol' are treated as zero.
%
%   See also SVD, IMA_SVD, PINV_SVD.

% F. Celi and F. Pasqualetti 2022

if ~exist('tol','var')
	tol = 1e-9;
end

[U S V] = svd(A);
r = min(size(S));

while r>0 && abs(S(r,r))<tol
    S(r,r) = 0;
    r = r - 1;
end

Si = inv(S(1:r,1:r));
S(1:r,1:r) = Si;
Z = V(:,r+1:end);