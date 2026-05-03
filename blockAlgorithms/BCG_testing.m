addpath("src")
addpath("matrix")

m = 24;
coef = 1;
reo = 1;
tol = 1e-13;
rankrevtol = 1e-12;
tillConv = 1;

%%% Random matrix %%%
% n = 500;
% A = rand(n,n);
% A = A + A';

%%% Matrix in /matrix %%%
load("bcsstk15.mat")
A = Problem.A;
[n,~] = size(A);

maxit = ceil(n/m)+20;

x_0 = zeros(n,m);
xex = [eye(m);zeros(n-m,m)];
b = A*xex;

tic
[x,T,omega,numrank,conv] = BCG(A,b,x_0,maxit,tillConv,coef,reo,tol,rankrevtol,xex);
fprintf("Done in %.4f s\n",toc);
fprintf('true convergence (time,iter) = (%f, %d)\n', conv(1), conv(2));
fprintf("fro || x - xex || = %.5e\n",norm(x - xex,"fro"))
fprintf("fro || A*x - b || = %.5e\n",norm(A*x - b,"fro"))


tic
[V,T_lan,numrank_lan] = BlockLanczos(A,b,maxit,coef,reo,rankrevtol);
fprintf("Lanczos done in %.4f s\n",toc);

T_lan = T_lan(1:size(T,1),1:size(T,2)); % in case of breakdown
fprintf("fro || T_BCG - T_Lan || = %.5e\n",norm(T_lan - T, "fro"))
diffs = computeDiffs(T_lan,T,m);


%% plot
% Block tridiagonal matrix
figure, hold on
plot(diffs)

title("Block tridiagonal matrices")
ylabel("$\Delta_k$", Interpreter="latex")
xlabel("Iterations")
set(gca, 'YScale', 'log')
set(gca, 'YMinorTick', 'off')

grid off
theme light
hold off

% Convergence
figure, hold on
plot(omega)

title("Convergence characteristics")
ylabel("$\omega_k$")
xlabel("Iteration")
set(gca, 'YScale', 'log')
set(gca, 'YMinorTick', 'off')

grid off
theme light
hold off

% Numerical ranks
figure, hold on
plot(numrank_lan(:,2),numrank_lan(:,1))
plot(numrank(:,2),numrank(:,1))

title("Numerical ranks")
ylabel("$rank_{\varepsilon}$",Rotation=0,Interpreter="latex")
legend("Block Lanczos","BCG")

grid off
theme light
hold off
