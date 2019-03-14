% Runs case study with gpr on bleen to calculate bootstrap estimates
clear;
clc;
%control panel
load dftips.mat
fid = fopen('log.txt', 'a');

level = 0.05; %confidence level

dftips = dftips(dftips.pkp_day.Month==1, :); % keep january data only
xIn = dftips.fare-15;yIn = dftips.tip_frac;% recenter around theshold
xMin = min(xIn); xMax= max(xIn);
nboot = 100;
nSamples =10000;
subsampleIndsL = nan(nSamples, nboot);
subsampleIndsR = nan(nSamples, nboot);

hyperoptstruct = struct('AcquisitionFunctionName','expected-improvement-plus');
hyperoptstruct.showPlots =false;
hyperoptstruct.Verbose =0;
hyperoptstruct.MaxObjectiveEvaluations =20;
hyperoptstruct.UseParallel = true;

gpr =struct();
[gprsLeft, gprsRight] = deal(cell(nboot, 1));

% pool = parpool(2);
%

load('case_study_gpr.mat') ;
for iB= 74:nboot
    
    fprintf('%d\n', iB)
    delete(gcp('nocreate'));
    %fprintf(fid, '%d --  %s\n', iB, datetime('now'));

    x = xIn; y = yIn;
    xL = x(x<0);yL = y(x<0);
    xR = x(x>=0);yR = y(x>=0);
    subsampleIndsL(:, iB) = randsample(length(xL), nSamples);
    subsampleIndsR(:, iB) = randsample(length(xR), nSamples);
    xL = xL(subsampleIndsL(:, iB)); yL= yL(subsampleIndsL(:, iB));
    xR = xR(subsampleIndsR(:, iB)); yR= yR(subsampleIndsR(:, iB));
    
%     tic;
%     gprMdll = fitrgp(xL, yL,  'KernelFunction', 'squaredexponential', ...
%     'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%     struct('AcquisitionFunctionName','expected-improvement-plus'));
%     fprintf('Left gpr, time elapsed %.3f\n', toc);
    
    tic;
    gprMdll = fitrgp(xL, yL,  'KernelFunction', 'squaredexponential', ...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    hyperoptstruct);%
    fprintf('Left gpr, time elapsed %.3f\n', toc);
    
    
    gpr(iB).lsigmaL= gprMdll.KernelInformation.KernelParameters(1);
    gpr(iB).lsigmaF= gprMdll.KernelInformation.KernelParameters(2);

    tic;
%     gprMdlr = fitrgp(xR, yR,  'KernelFunction', 'squaredexponential',  ...
%     'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%     struct('AcquisitionFunctionName','expected-improvement-plus'));
%     fprintf('Right gpr, time elapsed %.3f\n', toc)   
    gprMdlr = fitrgp(xR, yR,  'KernelFunction', 'squaredexponential',  ...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    hyperoptstruct);
    fprintf('Right gpr, time elapsed %.3f\n', toc)   
%     
    gpr(iB).rsigmaL= gprMdlr.KernelInformation.KernelParameters(1);
    gpr(iB).rsigmaF= gprMdlr.KernelInformation.KernelParameters(2);
    
    
    gprsLeft{iB} =  gprMdll; gprsRight{iB} = gprMdlr;

    % LATE estimate
    gpr(iB).tau =   predict(gprMdlr, 0)- predict(gprMdll, 0);
    

    % optimal threshold
    fun1 = @(x) (predict(gprMdlr, x)-predict(gprMdll, x));
    [cOptGPR, ~, exitflag]  = fzero(fun1, 0);
    if exitflag~=1
        cOptGPR = nan;
    elseif cOptGPR>xMax
        cOptGPR=xMax;% keep inside bw
    elseif cOptGPR<xMin
        cOptGPR=xMin;
    end

     gpr(iB).cOptGPR = cOptGPR;
    fprintf('%d ----cOptGPR---:%.3f\n', iB, cOptGPR);    
%    fprintf(fid, '%d ----cOptGPR---:%.3f\n', iB, cOptGPR);
save('case_study_gpr.mat', 'gpr')
end   
   %% 
% 
close all;figure(22); ah = gca; co=ah.ColorOrder;
legH =[]; % legend handles
% plot data

[xmeanR, ymeanR, ystd, nPoints]  = meanbin(xR, yR);
l=scatter(xmeanR, ymeanR, round(nPoints), 'Marker', '.', 'DisplayName', 'data');legH = [legH l];
% 
hold on%%
[xmeanL, ymeanL, ystd, nPoints]  = meanbin(xL, yL);
l=scatter(xmeanL, ymeanL, round(nPoints/1), 'Marker', '.', 'DisplayName', 'data');legH = [legH l];
% %
l =line([0, 0], ah.YLim, 'linewidth', 1, 'color', 'k', 'DisplayName','c');legH = [legH l];
ylim = ah.YLim;ah.XTickLabel= {'5','10','15','20','25'};

% plot gpr
% plot  extrapolation
xPredL = linspace(-10, 0)';
xPredR = linspace(0, 10)';
[yPredL,sL] = predict(gprMdll, xPredL);
[yPredR,sR] = predict(gprMdlr, xPredR);
[yPredRE, sRE] = predict(gprMdlr, xPredL);
[yPredLE, sLE] = predict(gprMdll, xPredR);


l=plot(xPredL, yPredL, 'linewidth', 1.5,  'DisplayName','GPR', 'color', co(4, :));legH = [legH l];gprColor= l.Color;
plot(xPredR, yPredR,'linewidth',1.5, 'color', gprColor);
plot(xPredR, yPredR+sR, 'color', gprColor);
plot(xPredL, yPredL+sL,'color', gprColor);
plot(xPredR, yPredR-sR, 'color', gprColor);
plot(xPredL, yPredL-sL,'color', gprColor);

plot(xPredL, yPredRE, 'color' ,gprColor, 'linestyle', '-.');
plot(xPredR, yPredLE, 'color' ,gprColor, 'linestyle', '-.');

plot(xPredL, yPredRE+sRE, 'color', gprColor);
plot(xPredL, yPredRE-sRE, 'color', gprColor);
plot(xPredR, yPredLE+sLE, 'color', gprColor);
plot(xPredR, yPredLE-sLE, 'color', gprColor);
% % pause;
