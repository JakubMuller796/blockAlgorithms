function [] = exp2_compar_plotting(matrix,ems,lw,gtol,plotdef)

fig = figure('Units','centimeters','Position',[2 2 22 16]);
tiledlayout(2,2,'Padding','loose','TileSpacing','loose');

for h = 1:4
m = ems(h);
filename = sprintf('exp/exp2/exp2_' + matrix + '_%d.mat',m);
load(filename);
den = real(trace(xex'*A*xex));

nexttile
hold on

% DRBCG
y_dr = computeNorms(omega_dr,den);
y_dr(y_dr < gtol) = gtol;
plot(y_dr,Color="red",LineWidth = lw, LineStyle="--")

% BCG
y_og = computeNorms(omega_bcg,den);
y_og(y_og < gtol) = gtol;
plot(y_og,Color="black",LineWidth = lw)

% BFBCG
y_bf = computeNorms(omega_bf,den);
y_bf(y_bf < gtol) = gtol;
plot(y_bf, Color = 'blue', LineWidth = lw,LineStyle="-.")

if size(numrank_bf,1) > 1 && plotdef
    change = numrank_bf(2:end,2);
    plot(change, y_bf(change), 'o',Color="blue", MarkerFaceColor='none',LineWidth=1,MarkerSize=4)

    yl = ylim;
    plot([numrank_bf(2,2) numrank_bf(2,2)], [y_bf(numrank_bf(2,2)) yl(1)], ...
         ':', LineWidth=0.4,Color="black");

elseif h == 1 && plotdef
    plot(nan, nan, 'o', ...
    'MarkerFaceColor','none', ...
    Color="blue", ...
    LineWidth=1);
end


ylabel('$\omega_k$',FontSize=12,Rotation=0,Interpreter="latex")
ylim([gtol*10^(-1) 1e+0])

xlabel('Počet iterací',FontSize=10)
xlim([0 maxits(h)])

ax = gca;
ax.XAxis.Exponent = 0;
xtickformat('%.0f')
xl = xlim;
xticks(linspace(xl(1), xl(2), 4))

%%% Uncomment to add true convergence iterations on the x axis.
%%% Unlike the graphs, in algorithms, we check for
%%% norm(r,"fro")/norm(b,"fro).
%%% So the iterations might be slightly off.

% newticks = sort([xticks, conv_bcg(2),conv_dr(2),conv_bf(2)]);
% newticks = unique(newticks);
% xticks(newticks)

if h == 1
    legend("DRBCG","BCG","BFBCG","deflation",Location="southwest")
end

title(sprintf("$m = %d$",ems(h)),FontWeight="normal",Interpreter="latex",FontSize=12)
set(gca, 'YScale', 'log')
set(gca, 'YMinorTick', 'off')

grid off
theme light
hold off

end