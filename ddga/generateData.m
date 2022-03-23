function data = generateData(sys,T,N,method,noise)
%GENERATEDATA Generate data experiments for a given discrete LTI system.
%   data = GENERATEDATA(sys,T,N,method,noise)
%   'method' :
%      - 'multiple': generate N experiments of length T
%      - 'hankel' : generate a single trajectory of length T
%   sys = (A,B,C) is the system model
%   If N is not specified N = n + m*T
%
%   data : ('method' = 'multiple')
%      - X0 initial state 
%      - X  full state trajectory
%      - U  inputs
%      - Y  outputs
%
%   data : ('method' = 'hankel')
%      - X  sates [x(0) x(1) ... x(T-1)]
%      - X1 sates [x(1) x(2) ... x(T)  ]
%      - U  inputs
%      - Y  outputs

if exist('noise','var')
    sigmaY = noise.sigmaY;
	sigmaX = noise.sigmaX;
	sigmaU = noise.sigmaU;
else
	sigmaY = 0;
	sigmaX = 0;
	sigmaU = 0;
end

% if exist('ic','var')
% 	sigmaX0 = ic;
% else 
% 	sigmaX0 = 10;
% end

O_x = observability_X(sys,T);
O_y = observability(sys,T);
F_x = forcing_X(sys,T);
F_y = forcing(sys,T);

n = size(sys.A,1);
m = size(sys.B,2);

data.method = method;
data.n = n;
data.m = m;
data.T = T;

switch method
	case 'multiple'
		% Number of experiments
		if size(N) == 0
			N = n+m*T;
		end

		X0 = normrnd(0,10,[n,N]);
		U  = normrnd(0,10,[m*T,N]);
		X = O_x * X0 + F_x * U;
		X = [X0; X];
		Y = O_y * X0 + F_y * U;
		
		data.X0noiseless = X0;
		data.Xnoiseless = X;
		data.Unoiseless = U;
		data.Ynoiseless = Y;

		data.X = X + normrnd(0, sigmaX, size(X));
		data.X0 = data.X(1:n,:);
		data.U = U + normrnd(0, sigmaU, size(U));
		data.Y = Y + normrnd(0, sigmaY, size(Y));

	case 'hankel'
		X0 = randn(n,1);
		U  = randn(m*T,1);

		X = O_x * X0 + F_x * U;
		X = [X0; X];
		Y = O_y * X0 + F_y * U;
		
		u = reshape(U,m,[]);
		U01T = u(:,1:end);
		x = reshape(X,n,[]);
		X0T = x(:,1:end-1);
		X1T = x(:,2:end);

		data.U = U01T;
		data.X = X0T;
		data.X1 = X1T;
		data.Y = Y;

	otherwise 
		println('Unknown option');
end