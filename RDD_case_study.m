clear;
clc;
addpath('util')
%control panel
load dftips.mat
level = 0.05;

dftips = dftips(dftips.pkp_day.Month==1, :); % keep january data only
x = dftips.fare-15;y = dftips.tip_frac;% recenter around theshold
%x = x(y<200); y = y(y<200);

%% LLR with uniform kernel in IK bandwidth
[est,se,h_opt]=rd_optbandwidth_uni(y, x, x>0,[],0,false,0);

indsIK = x>-h_opt & x<h_opt; % keep only values within bandwidht
xIK = x(indsIK); yIK = y(indsIK);
xIKL = xIK(xIK<0); xIKR = xIK(xIK>0); % split left and right data set
yIKL = yIK(xIK<0); yIKR = yIK(xIK>0); 


[llrMdl, llrMdlStats] = polyfit(xIKL, yIKL, 1);
[llrMdr, llrMdrStats] = polyfit(xIKR, yIKR, 1);


confIntL = conf_int_llr(linspace(-h_opt, h_opt), xIKL, yIKL, llrMdl, level);
confIntR =  conf_int_llr(linspace(-h_opt, h_opt), xIKR, yIKR, llrMdr, level);

%% GP with original bandwidth

xL = dftips(dftips.fare<15, :);
xR = dftips(dftips.fare>15, :);

tic;
gprMdlr_opt = fitrgp(xR(randsample(height(xR), 1000), 1:2), 'tip_frac','KernelFunction','squaredexponential',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));
gprMdlr = fitrgp(xR(:, 1:2), 'tip_frac',  'KernelFunction', 'squaredexponential', 'KernelParameters', gprMdlr_opt.KernelInformation.KernelParameters);
t =toc;
fprintf('Right gpr, time elapsed %.3f\n', toc)


tic;
gprMdll_opt = fitrgp(xL(randsample(height(xL), 1000), 1:2), 'tip_frac','KernelFunction','squaredexponential',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));
gprMdll = fitrgp(xL(:, 1:2), 'tip_frac',  'KernelFunction', 'squaredexponential', 'KernelParameters', gprMdll_opt.KernelInformation.KernelParameters);
t =toc;
fprintf('Left gpr, time elapsed %.3f\n', toc)

%% find optimal thresholds

load gprmodels.mat;
fun  = @(x) pi_x(x, gprMdll, gprMdlr);
cGPR = fzero(fun, 15);

tmp =llrMdr-llrMdl;
cLLR =- tmp(2)/tmp(1);


fun = @(x) tau_low(x, llrMdl, llrMdr, xIKL, xIKR, yIKL, yIKR, level);
cLLRCons = fzero(fun, 0);
[~, se, tau] =tau_low(0, llrMdl, llrMdr, xL, xR, yL, yR, level);
sprintf('%.3f & %.3f & %.3f  & %.3f  & %.3f & %.3f & %.3f & %.3f', llrMdl(1), llrMdl(2), llrMdr(1), llrMdr(2), h_opt, tau, cLLR+15, cLLRCons+15)

%% plot llr
clear; load trash/casestudy.mat;
close all;

figure(11); ah = gca;  co=ah.ColorOrder;
legH =[]; % legend handles
% plot data

[xmean, ymean, ystd, nPoints]  = meanbin(x, y);
l=scatter(xmean, ymean, round(nPoints/50), 'Marker', '.', 'DisplayName', 'data');legH = [legH l];
% 
hold on%%

l =line([0, 0], ah.YLim, 'linewidth', 1, 'color', 'k', 'DisplayName','c');legH = [legH l];
ylim = ah.YLim;
% plot KI optimal bandwidth
l =line([h_opt, h_opt], ah.YLim,'linestyle',  '--', 'linewidth', 1, 'color', 'k', 'DisplayName','IK');legH = [legH l];
line([-h_opt, -h_opt], ah.YLim, 'linestyle', '--', 'linewidth', 1, 'color', 'k');
% plot LLR extrapolation
l=plot([-h_opt, 0], polyval(llrMdl,[-h_opt, 0]), 'linewidth', 1.5, 'DisplayName','LLR', 'color', co(4, :));legH = [legH l];llrColor =l.Color;
plot([0, h_opt], polyval(llrMdr,[0, 0+h_opt]),'linewidth',1.5, 'color', llrColor);
plot([0-h_opt, 0], polyval(llrMdr,[0-h_opt, 0]), 'linewidth', 1.5, 'color' ,llrColor, 'linestyle', '-.');
%plot([0, h_opt], polyval(llrMdl,[0, h_opt]), 'linewidth', 1.5, 'color' ,llrColor, 'linestyle', '-.');

hf = fill_between(linspace(-h_opt, 0, 50), confIntL(1:50, 1), confIntL(1:50, 2));hf.FaceColor = [.9 .9 .9];
hf = fill_between(linspace(-h_opt, h_opt), confIntR(:, 1), confIntR(:, 2));hf.FaceColor = [.9 .9 .9];
ah.YLim=ylim;ah.XTickLabel= {'5','10','15','20','25'};


% plot optimal thresholds

line([cLLR, cLLR], ah.YLim,'linestyle',  '-.', 'linewidth', 1, 'color', llrColor, 'DisplayName', 'c^\star');%*_{LLR}');
line([cLLRCons, cLLRCons], ah.YLim, 'linestyle',  ':', 'linewidth', 1, 'color', llrColor, 'DisplayName', 'c^\star_{0.05}');%_{LLR, 0.05}');

l =scatter(cLLR, polyval(llrMdl, cLLR), 50, 'd',  'DisplayName','c*', 'MarkerFaceColor' , llrColor, 'MarkerEdgeColor' , llrColor);legH = [legH l];
l =scatter(cLLRCons, polyval(llrMdl, cLLRCons), 50, 's',  'DisplayName','c*_{0.05}', 'MarkerFaceColor' , llrColor, 'MarkerEdgeColor' , llrColor);legH = [legH l];

leg =legend(legH);
leg.FontSize=14;
xlabel(ah,'Fare($)', 'FontSize', 18)
ylabel(ah,'Tip percentage', 'FontSize', 18)

saveas(gcf, 'case_study_llr.png')
%% plot GPR models

figure(22);
xtest = 5:0.1:25;
xtest_r = 15:0.1:25;
xtest_l = 5:0.1:15;


l =plot(xtest_l-15,  predict(gprMdll, xtest_l'), 'linewidth', 1.5, 'color' , cm(3, :), 'DisplayName','GPR');legH = [legH l];
plot(xtest_l-15,  predict(gprMdlr, xtest_l'), 'linewidth', 1.5, 'color' , cm(3, :), 'linestyle', '-.');
plot(xtest_r-15,  predict(gprMdlr, xtest_r'), 'linewidth', 1.5, 'color' , cm(3, :));

% plot optimal thresholds
%l5 =scatter(cGPR, predict(gprMdll, cGPR-15), 50, 'd',  'DisplayName','c^*_{GPR}', 'MarkerFaceColor' , cm(3, :), 'MarkerEdgeColor' , cm(3, :));
l =scatter(cLLR, polyval(llrMdl, cLLR), 50, 'd',  'DisplayName','c^*_{LLR}', 'MarkerFaceColor' , cm(1, :), 'MarkerEdgeColor' , cm(1, :));legH = [legH l];
l =scatter(cLLRCons, polyval(llrMdl, cLLRCons), 50, 'd',  'DisplayName','c^*_{LLR}', 'MarkerFaceColor' , cm(4, :), 'MarkerEdgeColor' , cm(4, :));legH = [legH l];
% 
leg =legend(legH);
leg.FontSize=14;
xlabel(ah,'Fare($)', 'FontSize', 18)
ylabel(ah,'Tip percentage', 'FontSize', 18)

saveas(gcf, 'case_study.png')
