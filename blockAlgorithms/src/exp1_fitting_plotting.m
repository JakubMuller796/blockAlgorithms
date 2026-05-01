function [] = exp1_fitting_plotting(lw,styles,matrix,brk)

fig = figure('Units','centimeters','Position',[2 2 22 16]);
tiledlayout(2,1,'Padding','loose','TileSpacing','loose');

for mat = 1:2
distext = 0;

nexttile
hold on

for h = 1:3
filename = sprintf('exp/exp1/exp1_' + matrix(mat) + '_reo' + (h-1) + '.mat');
load(filename);

plot(diffs_dr,color="red",LineStyle=styles(h),LineWidth=lw)

plot(diffs_bcg,color="black",LineStyle=styles(h),LineWidth=lw)

if length(diffs_bcg) < length(diffs_dr) && brk
    x_last = length(diffs_bcg);
    y_last = diffs_bcg(end);

    plot([x_last x_last], [y_last (diffs_bcg(1)*1e-1)], ...
         ':', LineWidth=0.4,Color="black");
    
    if abs(distext - x_last) < 3
        text(x_last, diffs_bcg(1), num2str(x_last), ...
        'Rotation', 0, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontSize', 8);
    else
        text(x_last, diffs_bcg(1)*1e-1, num2str(x_last), ...
        'Rotation', 0, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontSize', 8);
        distext = x_last;
    end
end


set(gca, 'YScale', 'log')
xlabel("Počet iterací")
ylabel("$\Delta_k$",Interpreter="latex",FontSize=10,Rotation=0)
ylim([diffs_bcg(1)*1e-1, 10])
xlim([1 ceil(n/m)])

xt = xticks;
if ismember(ceil(n/m),xt) == false
    x_last = xt(end);
    if abs(x_last - ceil(n/m)) < abs(xt(1) - xt(2))/2
        xt = xt(1:end-1);
        xt = [xt, ceil(n/m)];
        xticks(sort(xt));
    else
        xt = [xt, ceil(n/m)];
        xticks(sort(xt));
    end
end

end

title("$\mathtt{" + matrix(mat) + "}$",FontWeight="normal",FontSize=12,Interpreter="latex")
if mat == 2 % legenda do druheho grafu
    legend("DRBCG (reo = 0)","BCG (reo = 0)","DRBCG (reo = 1)","BCG (reo = 1)", ...
        "DRBCG (reo = 2)", ...
        "BCG (reo = 2)",Location="southeast")
end

grid off
theme light
hold off

end