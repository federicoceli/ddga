function [A,B,C,sys] = randomLTI(n,m,p,stable)
% RANDOMLTI		Generates a random Linear Time-Invariant system.
%   [A,B,C,sys] = RANDOMLTI(n,m,p,stable)
%   n is the size of the system, m is the number of inputs and p is outputs
%   By deafult the system is Schur stable.

% F. Celi and F. Pasqualetti 2022

if ~exist('stable','var')
	stable = 1;
end

sys = drss(n,p,m);

A = sys.A;
B = sys.B;
C = sys.C;

sys = ss(A,B,C,[],-1);

if stable == 0
	A = A*2;
	sys.A = A;
end


end