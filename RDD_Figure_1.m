clear; close all;
addpath('util')
h =10; % bandwidth
c =0; % threshold
h_1 = c-h;
h_2 = c+h;

% left coefficients
a_0 = -10;
b_0 = 1.2;

fh = figure('PaperSize', [16, 12]);ah =gca;grey =[0.5, 0.5, 0.5];
ab_s  =[5, 3;  1.2 10;  -0.3, 12]    ;
colors = {'b', 'g', 'r', 'c'};
titles = {'\beta_1>\beta_0', '\beta_1=\beta_0', '\beta_1<\beta_0'};

for iter = 1:size(ab_s, 1)
    % plot outcome vs running variable
    subplot(3,3, iter); 
    ax= gca;hold on;ax.YLim=[-30, 55];ax.XLim =[-10 10];
    %// remove outer border
    box off %
    hold on
    a = axis; %// get axis size
    plot([a(1) a(2)],[a(3) a(3)],'k');
    plot([a(1) a(1)],[a(3) a(4)],'w', 'linewidth', 1.2); 
    % // end remove outer border
    plot([h_1,c], [a_0+b_0*h_1,a_0+b_0*c], 'linewidth',2, 'color', 'k');
    plot([c,h_2], [a_0+b_0*c,a_0+b_0*h_2], 'linewidth', 2, 'color','k', 'linestyle', ':');
 
    plot([c, c], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', ':');
    plot([h_1, h_1], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', '--');
    plot([h_2, h_2], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', '--');
    b_1= ab_s(iter, 1); a_1 = ab_s(iter, 2);
    plot([c,h_2], [a_1+b_1*c,a_1+b_1*h_2], 'linewidth',1.5, 'color',colors{iter});
    plot([h_1,c], [a_1+b_1*h_1,a_1+b_1*c], 'linewidth',1.5, 'color',colors{iter}, 'linestyle' , ':');
    ax.XTick =[-10,0, 10];
    ax.XTickLabel =[];%(['', '', ''])
    ax.YTick=[0];
    ax.YTickLabels = ['0'];
    ax.YAxis.FontSize =15;
    title(titles{iter}, 'FontSize', 20, 'color', colors{iter});
    if mod(iter, 3) ==1
        ylabel('Outcome', 'FontSize', 18);    
    end

    % plot treatment effect vs running variable
    subplot(3,3, iter+3); ax= gca;hold on;ax.YLim=[-30, 55];ax.XLim =[-10 10];
    %// remove outer border
    box off %
    hold on
    a = axis; %// get axis size
    plot([a(1) a(2)],[a(3) a(3)],'k');
    plot([a(1) a(1)],[a(3) a(4)],'w', 'linewidth', 1.2); 
    % // end remove outer border
    a_t=a_1-a_0;
    b_t=b_1-b_0;
    
    plot([c, c], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', ':');
    plot([h_1, h_1], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', '--');
    plot([h_2, h_2], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', '--');
    plot(ax.XLim, [0 0], 'linewidth', 1, 'color', grey, 'linestyle', '--');
    plot([h_1,h_2], [a_t+b_t*h_1,a_t+b_t*h_2], 'linewidth', 1.5, 'color', colors{iter});
    ar =fill_between([h_1, h_2], [0,0], [a_t+b_t*h_1,a_t+b_t*h_2]); ar.FaceColor = colors{iter};
    ax.XTick =[-10,0, 10];
    ax.XTickLabel =[];%(['', '', ''])
    ax.YTick=[0];
    ax.YTickLabels = ['0'];
    ax.YAxis.FontSize =15;
    if mod(iter, 3) ==1
        ylabel('Treatment effect', 'FontSize', 18);    
    end
    lgd = legend(ar,'\int_{c}^{c+h}{\pi(x)p(x)dx}', 'interpreter', 'latex');lgd.FontSize =8;

    
   % plot utility vs running variable
    subplot(3,3, iter+6); ax= gca;hold on;ax.YLim=[-120, 500];ax.XLim =[-10 10];
    %// remove outer border
    box off %
    hold on
    a = axis; %// get axis size
    plot([a(1) a(2)],[a(3) a(3)],'k');
    plot([a(1) a(1)],[a(3) a(4)],'w', 'linewidth', 1.2); 
    % // end remove outer border
    a_t=a_1-a_0;
    b_t=b_1-b_0;
    
    plot([c, c], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', ':');
    plot([h_1, h_1], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', '--');
    plot([h_2, h_2], ax.YLim, 'linewidth', 1, 'color', grey, 'linestyle', '--');
    plot(ax.XLim, [0 0], 'linewidth', 1, 'color', grey, 'linestyle', '--');%     x=np.arange(h_1, h_2, 0.01)
    x =h_1:0.01:h_2;
    
    plot(x,y_int(x, a_t, b_t, h_2), colors{iter})
%    
    xlabel('running variable', 'fontsize', 20);
    if iter==1
        scatter((-a_t/b_t), y_int(-a_t/b_t,  a_t, b_t, h_2), 40, 'MarkerFaceColor', colors{iter}, 'marker' ,'p')
        t = text((-a_t/b_t)+0.1, y_int(-a_t/b_t,  a_t, b_t, h_2)+20, 'c*', 'color', colors{iter});
    else
       scatter(h_1, y_int(h_1,  a_t, b_t, h_2), 40, 'MarkerFaceColor', colors{iter}, 'marker' ,'p')
       t = text(h_1+0.1, y_int(h_1,  a_t, b_t, h_2)+20, 'c*', 'color', colors{iter});    
    end
        
end
%saveas(fh, 'C:/Users/sofia/Dropbox/RDD/linear_cases', 'png')
