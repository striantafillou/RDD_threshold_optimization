% Runs case study with gpr (on bleen) to calculate bootstrap estimates

clear;
clc;
%control panel
load dftips.mat
level = 0.05; %confidence intervals

dftips = dftips(dftips.pkp_day.Month==1, :); % keep january data only
xIn = dftips.fare-15;yIn = dftips.tip_frac;% recenter around theshold
nboot = 100;
[~, bootInds] = bootstrp(nboot, @foo, xIn);


llr =struct();

for iB =1:nboot
    fprintf('%d\n', iB)
    x = xIn(bootInds(:, iB)); y = yIn(bootInds(:, iB));
    % LLR with uniform kernel in IK bandwidth
    [est,se,h_opt]=rd_optbandwidth_uni(y, x, x>0,[],0,false,0);
    indsIK = x>-h_opt & x<h_opt; % keep only values within bandwidth
    xIK = x(indsIK); yIK = y(indsIK);
    xIKL = xIK(xIK<0); xIKR = xIK(xIK>0); % split left and right data set
    yIKL = yIK(xIK<0); yIKR = yIK(xIK>0); 


    [llrMdl, llrMdlStats] = polyfit(xIKL, yIKL, 1);
    [llrMdr, llrMdrStats] = polyfit(xIKR, yIKR, 1);
    [~, se, tau] =tau_low(0, llrMdl, llrMdr, xIKL, xIKR, yIKL, yIKR, level);
    confIntL = conf_int_llr(linspace(-h_opt, h_opt), xIKL, yIKL, llrMdl, level);
    confIntR =  conf_int_llr(linspace(-h_opt, h_opt), xIKR, yIKR, llrMdr, level);

    tmp =llrMdr-llrMdl; % find optimal thresholds
    cLLR =- tmp(2)/tmp(1);

    fun = @(x) tau_low(x, llrMdl, llrMdr, xIKL, xIKR, yIKL, yIKR, level);
    cLLRCons = fzero(fun, 0);

    llr(iB).h_opt = h_opt;
    llr(iB).cLLR = cLLR;
    llr(iB).cLLRCons = cLLRCons;
    llr(iB).se = se;
    llr(iB).confIntL = confIntL;
    llr(iB).confIntR = confIntR;
    llr(iB).tau = tau;
    llr(iB).modelL = llrMdl;
    llr(iB).modelR = llrMdr;
    if iB ==1
        llr(nboot) = llr(iB);
    end
    if mod(iB, 10)==1
        save casestudy_llr.mat
    end
end

save casestudy_llr.mat
%% Results summary
load casestudy_llr.mat
coefs_0 = reshape([llr.modelL], 2, nboot);mc0 = mean(coefs_0,2);std0 = std(coefs_0,[], 2);
coefs_1 = reshape([llr.modelR], 2, nboot);mc1 = mean(coefs_1,2);std1 = std(coefs_1,[], 2);
sprintf('%.3f & %.3f & %.3f  & %.3f  & %.3f & %.3f & %.3f & %.3f', mc0(1), mc0(2), mc1(1), mc1(2), mean([llr.h_opt]),mean([llr.tau]), mean([llr.cLLR])+15, mean([llr.cLLRCons])+15)
sprintf('(%.3f) & (%.3f) & (%.3f)  & (%.3f)  & (%.3f) & (%.3f) & (%.3f) & (%.3f)', std0(1), std0(2), std1(1), std1(2), std([llr.h_opt]),std([llr.tau]), std([llr.cLLR])+15, std([llr.cLLRCons])+15)
