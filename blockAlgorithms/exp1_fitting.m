addpath("src")
addpath("matrix")

matrix = ["bcsstk15","crystm01"];

%%% exp1 setup %%%
m = 24;
coef = 1;
tol = 1e-13;
rankrevtol = 1e-12;
tillConv = 0;
%%% exp1 setup %%%

for mat = 1:2 
fprintf("==== " + matrix(mat) + " ====\n")

load(matrix(mat) + ".mat")
A = Problem.A;
fprintf("condition number: %.e\n",condest(A))

n = size(A,1);
maxit = ceil(n/m);
xex = [eye(m,m);zeros(n-m,m)];
x_0 = zeros(n,m);
b = A*xex;

times = zeros(3,1);

for h = 1:3
reo = h-1;
fprintf("==== reo = " + reo + " ====\n")

%%% block Lanczos
tic
[V_lan,T_lan,numrank_lan] = BlockLanczos(A,b,maxit,coef,reo,rankrevtol);
times(1) = toc;
fprintf("Lanczos done in %.4f s\n",times(1));
if size(numrank_lan,1) > 1
    fprintf("Lanczos first rank-deficiency in: %.d\n",numrank_lan(2,2))
end

%%% BCG
tic
[x_bcg, T_bcg,~,numrank_bcg,conv_bcg] = BCG(A,b,x_0,maxit,tillConv,coef,reo,tol,rankrevtol,[]);
times(2) = toc;
fprintf("BCG done in %.4f s\n",times(2));
fprintf("fro || T_lan - T_bcg || = %.5e\n",norm(T_bcg-T_lan(1:size(T_bcg,1),1:size(T_bcg,2)),"fro"))

%%% DRBCG
tic
[x,T_dr,~,numrank_dr,conv_dr] = DRBCG(A,b,x_0,coef,reo,tol,rankrevtol,maxit,tillConv,"qr",[]);
times(3) = toc;
fprintf("DR done in %.4f s\n",times(3));
fprintf("fro || T_lan - T_dr || = %.5e\n",norm(T_dr-T_lan,"fro"))

%%% diffs
fprintf("=== Algs done, computing diffs ===\n")

diffs_dr = computeDiffs(T_lan,T_dr,m);
diffs_bcg = computeDiffs(T_lan,T_bcg,m);

fprintf("=== Diffs done ===\n")

filename = sprintf('exp/exp1/exp1_' + matrix(mat) + '_reo' + (h-1) + '.mat');
save(sprintf(filename), ...
    "n","m","T_lan","T_bcg","T_dr","diffs_dr","diffs_bcg","times", ...
    "maxit","conv_bcg","conv_dr","numrank_dr","numrank_bcg","numrank_lan");

end
end

%% plot
%%% styles are the three lineStyles for reo = 0,1,2
%%% lw stands for line width
%%% brk will plot the last iteration of BCG

styles = ["-","--",":"];
lw = 1.2;
brk = 0;

exp1_fitting_plotting(lw,styles,matrix,brk)
