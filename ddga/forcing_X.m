function F_x = forcing_X(sys,T)
% FORCING_X Generates the state forcing matrix of sys = (A,B,C).
%   F_x = FORCING_X(sys,T)
%        _                          _
%       |   B          ...   0     0 |
%       |   AB         ...   0     0 |
% F_x = |   A^2B       ...   0     0 |
%       |   ....       ...   ...   0 | 
%       |   CA^(T-1)B  ...   AB    B |
%        ¯                          ¯

% F. Celi and F. Pasqualetti 2022

A = sys.A;
B = sys.B;
C = sys.C;

n = size(A,1);
F_x = zeros(n*T,size(B,2)*T);

for r=1:T+1
	for c=1:r-1
    	F_x((r-1)*n+1-n:r*n-n,(c-1)*size(B,2)+1:c*size(B,2))=(A^(r-1-c))*B;
	end
end
