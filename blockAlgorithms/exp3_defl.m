addpath("src")
addpath("matrix")


matrix = "bcsstk15";
load(sprintf(matrix + ".mat"))
A = Problem.A;
[n,~] = size(A);

%%% exp3 setup %%%
m = 32;
maxit = 6000;
rankrevtols = [1e-12,1e-8];
facs = ["svd","svd"];
tol = 1e-13;
%%% exp3 setup %%%

x_0 = zeros(n,m);
xex = [eye(m);zeros(n-m,m)];
b = A*xex;

times = zeros(2,3);

for h = 1:2

rankrevtol = rankrevtols(h);
fac = facs(h);

fprintf("==== orth = " + fac + ", rankrevtol = " + rankrevtol + " ====\n")

%%% BCG
tic
[x_bcg,~,bcg_omega,numrank_bcg_cell{h},conv_bcg_cell{h}] = BCG(A,b,x_0,maxit,0,0,tol,rankrevtol,xex);
times(h,1) = toc;
fprintf("BCG fro || x - xex || = %.5e\n",norm(x_bcg - xex,"fro"))
fprintf("BCG fro || A*x - b || = %.5e\n",norm(A*x_bcg - b,"fro"))
fprintf("BCG done in %.4f s\n",toc);
fprintf('true convergence (time,iter) = (%f, %d)\n', conv_bcg_cell{h}(1), conv_bcg_cell{h}(2));

%%% DRBCG
tic
[x_dr,~,dr_omega,numrank_dr_cell{h},conv_dr_cell{h}] = DRBCG(A,b,x_0,0,0,tol,rankrevtol,maxit,"qr",xex);
times(h,2) = toc;
fprintf("DR fro || x - xex || = %.5e\n",norm(x_dr - xex,"fro"))
fprintf("DR fro || A*x - b || = %.5e\n",norm(A*x_dr - b,"fro"))
fprintf("DR done in %.4f s\n",toc);
fprintf('true convergence (time,iter) = (%f, %d)\n', conv_dr_cell{h}(1), conv_dr_cell{h}(2));

%%% BFBCG
tic
[x_bf,omega_bf_cell{h},numrank_bf_cell{h},conv_bf_cell{h}] = BFBCG(A,b,x_0,"svd",tol,rankrevtol,maxit,xex);
times(h,3) = toc;
fprintf("BF fro || x - xex || = %.5e\n",norm(x_bf - xex,"fro"))
fprintf("BF fro || A*x - b || = %.5e\n",norm(A*x_bf - b,"fro"))
fprintf("BF with rankrevtol = %d done in %.4f s\n",rankrevtol,toc);
fprintf('true convergence (time,iter) = (%f, %d)\n', conv_bf_cell{h}(1), conv_bf_cell{h}(2));

end

save(sprintf('exp/exp3/exp3.mat'), ...
    "dr_omega","bcg_omega",'numrank_dr_cell', ...
    "omega_bf_cell","numrank_bf_cell", ...
    "numrank_bcg_cell",'xex','A','m','maxit','tol',"times", ...
    "rankrevtols","conv_bcg_cell","conv_bf_cell","conv_dr_cell");

times

%% plot

%%% lw stands for line width
%%% gtol stands for graph tolerence (usually the same as tol)
%%% brk will plot the last iteration of BCG

lw = 1.2;
gtol = 1e-13;
brk = 1;

exp3_defl_plotting(lw,gtol,brk)
