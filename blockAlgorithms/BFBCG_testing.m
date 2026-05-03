addpath("src")
addpath("matrix")

m = 10;
tol = 1e-15;
rankrevtol = 1e-12;
tillConv = 1;
type = "svd";

%%% Random matrix %%%
% n = 500;
% A = rand(n,n);
% A = A + A';

%%% Matrix in /matrix %%%
load("bcsstk15.mat")
A = Problem.A;
[n,~] = size(A);

maxit = 10000;

x_0 = zeros(n,m);
xex = [eye(m);zeros(n-m,m)];
b = A*xex;

tic
[x,bf_omega,numrank,conv] = BFBCG(A,b,x_0,type,tol,rankrevtol,maxit,tillConv,xex);
fprintf("Done in %.4f s\n",toc);
fprintf('true convergence (time,iter) = (%f, %d)\n', conv(1), conv(2));
fprintf("fro || x - xex || = %.5e\n",norm(x - xex,"fro"))
fprintf("fro || A*x - b || = %.5e\n",norm(A*x - b,"fro"))


%% plot
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
plot(numrank(:,2),numrank(:,1))

title("Numerical rank")
ylabel("$rank_{\varepsilon}$",Rotation=0,Interpreter="latex")
legend("BFBCG")

grid off
theme light
hold off
