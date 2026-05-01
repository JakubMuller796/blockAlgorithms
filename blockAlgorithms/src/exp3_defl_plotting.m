function [] = exp3_defl_plotting(lw,gtol,brk)

% fig convergence

fig = figure('Units','centimeters','Position',[2 2 22 8]);
tiledlayout(1,1,'Padding','loose','TileSpacing','loose');

filename = sprintf('exp/exp3/exp3.mat');
load(filename);

den = real(trace(xex'*A*xex));

nexttile
hold on

% DRBCG
y_dr = computeNorms(dr_omega,den);
y_dr(y_dr < gtol) = gtol;
h1 = plot(y_dr,Color="red",LineWidth = lw, LineStyle="--");

% BCG
y_bcg = computeNorms(bcg_omega,den);
y_bcg = y_bcg(~isnan(y_bcg));
y_bcg(y_bcg < gtol) = gtol;
h2 = plot(y_bcg,Color="black",LineWidth = lw);

if length(y_bcg) < maxit && brk
    x_last = length(y_bcg);
    y_last = y_bcg(end);

    plot([x_last x_last], [y_last (gtol*10^(-2))], ...
         ':', LineWidth=0.4,Color="black");
    
    text(x_last, gtol*10^(-1), num2str(x_last), ...
    'Rotation', 0, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'FontSize', 8);
end


% BFBCG
for tols = 1:length(rankrevtols)
    bf_omega = omega_bf_cell{tols};
    numrank_bf = numrank_bf_cell{tols};

    y_bf = computeNorms(bf_omega,den);
    y_bf(y_bf < gtol) = gtol;
    h3 = plot(y_bf, Color = 'blue', LineWidth = lw,LineStyle="-.");

    exp = round(log10(double(rankrevtols(tols))));
    desc = sprintf('$\\varepsilon = 10^{%d}$',exp);
    text(length(y_bf), y_bf(end), desc, ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex');  
end



xlabel('Počet iterací',FontSize=10)
ylabel('$\omega_k$',FontSize=12,Rotation=0,Interpreter="latex")

ylim([gtol*10^(-1) 1e+0])
xlim([0 maxit])

ax = gca;
ax.XAxis.Exponent = 0;
xtickformat('%.0f')

legend([h1 h2 h3],"DRBCG","BCG","BFBCG",Location="northeast")

title("$ m = " + m + "$",FontWeight="normal",FontSize=12,Interpreter="latex")

set(gca, 'YScale', 'log')
set(gca, 'YMinorTick', 'off')

grid off
theme light
hold off


%% fig numrank

fig = figure('Units','centimeters','Position',[2 2 (length(rankrevtols))*11 8]);
tiledlayout(1,length(rankrevtols),'Padding','loose','TileSpacing','loose');

for h = 1:length(rankrevtols)
    
    nexttile
    hold on

    % DRBCG
    numrank_dr = numrank_dr_cell{h};
    numrank_dr = [numrank_dr;[numrank_dr(end,1) maxit]]; 
    plot(numrank_dr(:,2),numrank_dr(:,1),'b', 'LineWidth',lw,Color="red")

    %BCG
    numrank_bcg = numrank_bcg_cell{h};
    plot(numrank_bcg(:,2),numrank_bcg(:,1),'b', 'LineWidth',lw,Color="black")
    
    % BFBCG
    numrank_bf = numrank_bf_cell{h};
    numrank_bf = [numrank_bf;[numrank_bf(end,1) maxit]];    
    plot(numrank_bf(:,2),numrank_bf(:,1),'b', 'LineWidth',lw,Color="blue",LineStyle="-.")
    
    ax = gca;
    yt = ax.YTick;
    ticks = sort(unique([yt, m]));
    yticks(ticks)
    ylim([-1,numrank_bf_cell{1}(1,1) + 1])

    exp = round(log10(double(rankrevtols(h))));
    title("$\varepsilon = 10^{" + exp + "}$", ...
          Interpreter="latex")    

    xlabel('Počet iterací',FontSize=10)
    ylabel('$\mathrm{rank}_{\varepsilon}$',Interpreter="latex",FontSize=12,Rotation=0)
    
    if h == 1
        legend("DRBCG","BCG","BFBCG",Location="southwest")
    end

    grid off
    theme light
    hold off
end