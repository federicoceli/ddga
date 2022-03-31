%%
close all
clear all
addpath('./utilities/')
tol = 1e-9;
set_graphics;

%% System setup 
n = 2;
m = 1;
p = 1;

A = [0.5 0 ; 0 1];
C = [3 -5];
C1 = [3 4];
B = [1; 0];
Sigma1 = ss(A,B,C,[]);
Sigma2 = ss(A,B,C1,[]);
T = 3;
N = 5;

%% Model-based geometric control invariants
Vstar_mb_S1 = vstar(Sigma1);
Vstar_mb_S2 = vstar(Sigma2);

%% Data-driven geometric control invariants (Sigma1)
% We show this step by step as an exercise.
% The same results can be obtain by running:
%
% data = generateData(sys,T,[],'multiple');
% Vstar_dd = vstardd(data,tol);

O_x = observability_X(Sigma1,T);
O_y = observability(Sigma1,T);
F_x = forcing_X(Sigma1,T);
F_y = forcing(Sigma1,T);

% Generate data as multiple (N) experiments of length T
X0 = [1 0 0 0 0;
	  0 1 0 0 0;];
% inputs
U = [0 0 1 0 0;
	 0 0 0 1 0;
	 0 0 0 0 1];
% state trajectoriesX
X = O_x * X0 + F_x * U;
X = [X0; X];
% output trajectories
Y = O_y * X0 + F_y * U;

data.n = n;
data.m = m;
data.X0 = X0;
data.X = X;
data.U = U;
data.Y = Y;

% <-- Verify this matches Vstar_mb for Sigma1 -->
Vstar_dd_S1 = vstardd(data,tol);

%% System Id (MOESP) + Traditional Geometric Control
% The system is not in minimal form and as a result traditional
% input-output system identification algorithm, such as MOESP cannot
% succesfully identify a model of the system. 

N = 200;
u = randn(N,1);
y = dlsim(A,B,C,[],u);

k = 15;
[ss1,ssfun] = moesp(y,u,k);
[A_id,B_id,C_id]=ssfun(2);


%% (Sigma2)
O_x1 = observability_X(Sigma2,T);
O_y1 = observability(Sigma2,T);
F_x1 = forcing_X(Sigma2,T);
F_y1 = forcing(Sigma2,T);

% Generate data as multiple (N) experiments of length T
X01 = [1 0 0 0 0;
	  0 -1.25 0 0 0;];
% inputs
U = [0 0 1 0 0;
	 0 0 0 1 0;
	 0 0 0 0 1];
% state trajectoriesX
X1 = O_x1 * X01 + F_x1 * U;
X1 = [X01; X1];
% output trajectories
Y1 = O_y1 * X01 + F_y1 * U;

data1.n = n;
data1.m = m;
data1.X0 = X01;
data1.X = X1;
data1.U = U;
data1.Y = Y1;

% <-- Verify this matches Vstar_mb for Sigma2 -->
Vstar_dd_S2 = vstardd(data1,tol);