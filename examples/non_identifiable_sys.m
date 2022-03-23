%%
close all
clear all
addpath('../ddga/')
tol = 1e-9;
set_graphics;

%% System setup 
n = 2;
m = 1;
p = 1;

A = [1 2; 2 4];
C = [2 4];
B = [1; 2];
sys = ss(A,B,C,[]);

T = 3;
N = 5;

%% Model-based geometric control invariants
[Vstar_mb, F_mb] = vstar(A,B,C);
Sstar_mb = sstar(A,B,C);
Rstar_mb = ints(Vstar_mb,Sstar_mb);

%% Data-driven geometric control invariants
% We show this step by step as an exercise. The same results can be obtain
% by running the following lines of code:
%
% data = generateData(sys,T,[],'multiple');
% Vstar_dd = vstardd(data,tol);

O_x = observability_X(sys,T);
O_y = observability(sys,T);
F_x = forcing_X(sys,T);
F_y = forcing(sys,T);

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

% <-- Verify this matches Vstar_mb! -->
Vstar_dd = vstardd(data,tol);

%% System Id (MOESP) + Traditional Geometric Control
% The system is not in minimal form and as a result traditional
% input-output system identification algorithm, such as MOESP cannot
% succesfully identify a model of the system. 

N = 200;
u = randn(N,1);
y = dlsim(A,B,C,[],u);

k = 15;
[ss,ssfun] = moesp(y,u,k);
[A_id,B_id,C_id]=ssfun(2);

% <-- Verify this does NOT match Vstar_mb! -->
Vstar_id = vstar(A_id,B_id,C_id); 

