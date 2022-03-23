function Q = ima_svd(A,tol)
% IMA_SVD	Orthogonalization.
%   Q = IMA_SVD(A,tol) is an orthonormal basis for imA.
%   The computation is based on SVD and any
%   singular value less than tolerance 'tol' are treated as zero.
%
%   See also SVD, PINV_SVD, NULL_SVD.

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

Si=inv(S(1:r,1:r));
S(1:r,1:r) = Si;
Q = ortco(ortco(ima(U*S*V')));

if sum(size(Q)) == 0
	Q = zeros(size(U,1),1);
end