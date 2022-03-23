function O_x = observability_X(sys,T)
% OBSERVABILITY_X Generates the state free response matrix of sys = (A,B,C).
%   O_x = OBSERVABILITY_X(sys,T)
%        _       _
%       |   A     |
%       |   A^2   |
% O_x = |   ...   |
%       |   CA^T  |
%        ¯       ¯

% F. Celi and F. Pasqualetti 2022

A = sys.A;
B = sys.B;
C = sys.C;
	
O_x = zeros(size(A,1)*T,size(A,2));

for r = 1:T
	O_x(1+(r-1)*size(A,1):size(A,1)+(r-1)*size(A,1),:) = A^(r);
end

