function Sstar = sstardd(data,tol)
% SSTARDD	Computes an estimate of the smallest input-containing 
%   conditioned invariant of the system that generated data.
%   Sstar = SSTARDD(data,tol)
%   Use generateData(A,B,C,D,T,N,'multiple',noise) to generate N
%   experiments of lenght T. All operations based on singular value
%   decomposition are run with tolerance tol.
%
%   See also VSTARDD, GENERATEDATA, NULL_SVD, IMA_SVD.

% F. Celi and F. Pasqualetti 2022

if ~exist('tol','var')
	tol = 1e-9;
end

Y = data.Y;
X = data.X;
X0 = data.X0;

n = size(X0,1);
Xf = X(end-n+1:end,:);

K0 = null_svd(X0, tol);
beta_s = null_svd(Y*K0,tol);

% For robustness issues remove small singular values
[Us,Ss,Vs] = svd(Xf*K0*beta_s,'econ');
for i = 1 : size(Ss,1)
	if Ss(i,i) < tol
		Ss(i,i) = 0;
	end
end

Sstar = ima_svd(Us*Ss*Vs',tol);