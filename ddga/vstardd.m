function [Vstar, F] = vstardd(data,tol)
% VSTARDD	Computes an estimate of the maximum output-nulling
%   controlled invariant of the system that generated data.
%   [Vstar, F] = VSTARDD(data,tol)
%   Use generateData(A,B,C,D,T,N,'multiple',noise) to generate N
%   experiments of lenght T. All operations based on singular value
%   decomposition are run with tolerance tol.
%
%   See also SSTARDD, GENERATEDATA, NULL_SVD, IMA_SVD.

% F. Celi and F. Pasqualetti 2022

if ~exist('tol','var')
	tol = 1e-9;
end

Y = data.Y;
X0 = data.X0;
U = data.U;

% Compute Vstar with data
Vstar = ima_svd(X0*null_svd(Y,tol),tol);

% Compute a friend F

V = Vstar;
n = data.n;
m = data.m;
X = data.X(n+1:end,:);

X0T = reshape(X(1:end-n,1),n,[]);
X1T = reshape(X(n+1:end,1),n,[]);
U0T = reshape(U(m+1:end,1),m,[]);

gamma = - pinv((eye(n) - V*pinv(V))*X1T*null(X0T)) * ...
	(eye(n) - V*pinv(V))*X1T*pinv(X0T)*V*pinv(V);

F = U0T*(pinv(X0T) + null(X0T)*gamma);

end