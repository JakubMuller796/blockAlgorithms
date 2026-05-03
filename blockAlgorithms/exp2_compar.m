addpath("src")
addpath("matrix")


matrix = "bcsstk15";
load(sprintf(matrix + ".mat"))
A = Problem.A;
[n,~] = size(A);

times = zeros(4,3);

%%% exp2 setup %%%
ems = [1,4,16,24];
maxits = [24000,12000,9000,6000];
tol = 1e-13;
rankrevtol = 1e-12;
tillConv = 0;
%%% exp2 setup %%%

for h = 1:4

m = ems(h);
maxit = maxits(h);

x_0 = zeros(n,m);
xex = [eye(m);zeros(n-m,m)];
b = A*xex;

fprintf("==== m = %d ====\n",m)

%%% BCG
tic
[x_bcg,~,omega_bcg,numrank_bcg,conv_bcg] = BCG(A,b,x_0,maxit,tillConv,0,0,tol,0,xex);
times(h,1) = toc;
fprintf("BCG fro || x - xex || = %.5e\n",norm(x_bcg - xex,"fro"))
fprintf("BCG fro || A*x - b || = %.5e\n",norm(A*x_bcg - b,"fro"))
fprintf("BCG done in %.4f s\n",toc);
fprintf('true convergence (time,iter) = (%f, %d)\n', conv_bcg(1), conv_bcg(2));

%%% DRBCG
tic
[x_dr,~,omega_dr,numrank_dr,conv_dr] = DRBCG(A,b,x_0,0,0,tol,0,maxit,tillConv,"qr",xex);
times(h,2) = toc;
fprintf("DR fro || x - xex || = %.5e\n",norm(x_dr - xex,"fro"))
fprintf("DR fro || A*x - b || = %.5e\n",norm(A*x_dr - b,"fro"))
fprintf("DR done in %.4f s\n",toc);
fprintf('true convergence (time,iter) = (%f, %d)\n', conv_dr(1), conv_dr(2));

%%% BFBCG
tic
[x_bf,omega_bf,numrank_bf,conv_bf] = BFBCG(A,b,x_0,"svd",tol,rankrevtol,maxit,tillConv,xex);
times(h,3) = toc;
fprintf("BF fro || x - xex || = %.5e\n",norm(x_bf - xex,"fro"))
fprintf("BF fro || A*x - b || = %.5e\n",norm(A*x_bf - b,"fro"))
fprintf("BF done in %.4f s\n",toc);
fprintf('true convergence (time,iter) = (%f, %d)\n', conv_bf(1), conv_bf(2));


save(sprintf('exp/exp2/exp2_' + matrix + '_%d.mat',m), ...
    "omega_bf","omega_dr","omega_bcg", ...
    'numrank_bcg','numrank_dr','numrank_bf','xex','A', ...
    'n','m','maxit','maxits','tol',"times","conv_bf","conv_dr","conv_bcg");

end

times
fprintf("==== Computations done, plotting ==== \n")

%% plot

%%% lw stands for line width
%%% gtol stands for graph tolerence (usually the same as tol)
%%% plotdef stands for plotting deflation points (only in case of BFBCG)

ems = [1,4,16,24];
lw = 1.2;
gtol = 1e-13;
plotdef = 0;

exp2_compar_plotting(matrix,ems,lw,gtol,plotdef)
