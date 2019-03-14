clear; close all;
addpath('util')

h =10; % bandwidth
c =0; % threshold
h_1 = c-h;
h_2 = c+h;




% simulate data;
nSamples =100; % number of samples
sigma_e =10;% error variance
x = -10+20*rand(nSamples, 1);  % simulate x in [-10,10]
e = sigma_e*randn(nSamples, 1);  % noise ~ N(0,10)
const = ones(nSamples,1);

%split left and right
iL= x<0; iR = x>=0; 
eL = e(iL); eR= e(iR);
xL =x(iL); xR = x(iR);
constL = const(iL); constR = const(iR);

% simulate yL = b0*xL+a0+eL, yR = b1*xR+a1+eR
a_0 = -10;
b_0 = 1.2;
ab_s  =[5, 3] ;
b_1= ab_s(1, 1); a_1 = ab_s(1, 2);
y = nan(nSamples,1);
yL = [xL constL]*[b_0; a_0]+eL;
yR = [xR constR]*[b_1; a_1]+eR;
y(iL) = yL; y(iR)=yR;

% fit llr models
[llrMdl, llrMdlStats] = polyfit(xL, yL, 1);
[llrMdr, llrMdrStats] = polyfit(xR, yR, 1);


% Plot data and optimal threshold
figure(11);grey =[0.5, 0.5, 0.5];colors ={'b'};
ax= gca;hold on;ax.YLim=[-30, 55];ax.XLim =[-10 10]; 

%// remove outer border
box off %
hold on
a = axis; %// get axis size
plot([a(1) a(2)],[a(3) a(3)],'k');
plot([a(1) a(1)],[a(3) a(4)],'w', 'linewidth', 1.2); 
% // end remove outer border

%plot vertical/horizontal lines
lc =plot([c, c], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', ':');
plot([h_1, h_1], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', '--');
plot([h_2, h_2], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', '--');
%plot(ax.XLim, [0 0], 'linewidth', 2, 'color', grey, 'linestyle', '-');


% scatter points;
scatter(xL, yL, 100, '.', 'MarkerEdgeColor', 'k');
scatter(xR, yR, 100, '.', 'MarkerEdgeColor', colors{1});

% plot fitted lines
plot([c,h_2], [polyval(llrMdr, c),polyval(llrMdr, h_2)], 'linewidth',1.5, 'color',colors{1});
plot([h_1,c],  [polyval(llrMdr, h_1),polyval(llrMdr, c)], 'linewidth',1.5, 'color',colors{1}, 'linestyle' , ':');
plot([h_1,c],  [polyval(llrMdl, h_1),polyval(llrMdl, c)], 'linewidth',1.5, 'color','k');

% labels etc
ax.XTick =[-10,0, 10];
ax.XTickLabel =[];%(['', '', ''])
ax.YTick=[0];
ax.YTickLabels = ['0'];
ax.YAxis.FontSize =15;

ylabel('Outcome', 'FontSize', 18, 'interpreter', 'latex');    
%ylabel('$\hat{\pi}(x)$', 'FontSize', 18, 'interpreter', 'latex');    
xlabel('$x$', 'FontSize', 18, 'interpreter', 'latex')

%% Figure 2
figure(22);linecolors = ax.ColorOrder;
ax= gca;hold on;ax.YLim=[-30, 55];ax.XLim =[-10 10]; 
%// remove outer border
box off %
hold on
a = axis; %// get axis size
plot([a(1) a(2)],[a(3) a(3)],'k');
plot([a(1) a(1)],[a(3) a(4)],'w', 'linewidth', 1.2); 
% // end remove outer border

%plot vertical/horizontal lines
plot([c, c], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', ':');
plot([h_1, h_1], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', '--');
plot([h_2, h_2], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', '--');
plot(ax.XLim, [0 0], 'linewidth', 2, 'color', grey, 'linestyle', '-');

x_new = -10:0.5:10;
levels = [0.1, 0.05, 0.01];
for iLevel =1:3
    level= levels(iLevel);
    [tau_lows, ses, taus] =tau_low(x_new,  llrMdl, llrMdr, xL, xR,  yL, yR, level);
    [cOptLLR(iLevel), cOptLLRCons(iLevel), output] = c_opt_llr(x, y, level, false, 10);

    if iLevel==1
      lh2(1) =  plot(x_new, taus, 'color', 'k', 'linewidth', 2);
      scatter(cOptLLR(iLevel), 0, 70, 'h', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
      figure(11); lh1(1)=plot([cOptLLR(iLevel), cOptLLR(iLevel)], ax.YLim,  'color', 'k', 'linestyle', ':', 'linewidth', 2); figure(22)
    end
    lh2(iLevel+1)=plot(x_new, tau_lows,  'linestyle','-', 'color', linecolors(iLevel, :), 'linewidth', 2);
    scatter(cOptLLRCons(iLevel), 0, 70, 'h', 'MarkerFaceColor',  linecolors(iLevel, :), 'MarkerEdgeColor', 'none');
    figure(11);lh1(iLevel+1) =plot([cOptLLRCons(iLevel), cOptLLRCons(iLevel)], ax.YLim, 'color', linecolors(iLevel, :), 'linestyle', ':', 'linewidth', 2); figure(22)
end
lg =legend(lh2,{'$\hat{\pi}(x)$', '$\alpha =0.1$', '$\alpha =0.05$', '$\alpha =0.01$'}, 'location', 'northwest', 'FontSize', 18);
lg.Interpreter ='latex';

% labels etc
ax.XTick =[-10,0, 10];
ax.XTickLabel =[];%(['', '', ''])
ax.YTick=[0];
ax.YTickLabels = ['0'];
ax.YAxis.FontSize =15;

ylabel('$\hat{\pi}(x)$', 'FontSize', 18, 'interpreter', 'latex');    
xlabel('$x$', 'FontSize', 18, 'interpreter', 'latex')
saveas(gcf, 'C:\Users\sofia\Dropbox\RDD\llr_cons_2.png')

figure(11);lg1=legend([lc lh1],{'$c$','$c^*$', '$c^*_{0.1}$', '$c^*_{0.05}$', '$c^*_{0.01}$'}, 'location', 'northwest', 'FontSize', 18);
lg1.Interpreter ='latex';
%saveas(gcf, 'C:\Users\sofia\Dropbox\RDD\llr_cons_1.png')
