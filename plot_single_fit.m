% function plot_data = plot_single_fit(subjidx,mflags)
% 
% Produces a plot with the maximum likelihood fit of the model specified by
% mflags to the data of the subject specified by subjidx. The flags of the 
% main model used in the paper (P-C-VP) are [3 3 2].
%
% INPUT
%  subjidx : subject index (integer in  range 1-6)
%  mflags  : 1x3 vector with model flags:
%            mflags(1): 1=All, 2=Fixed, 3=Poisson number of items encoded
%            mflags(2): 1=Quantized-even, 2=Quantized-random, 3=Continuous
%            mflags(3): 1=EP, 2=VP
%
% EXAMPLE
%  >> plot_single_fit(1,[3 3 2])  % Plot fit of main model (P-C-VP) to data of subject 1
%  
% This file is part of the code published with the paper "Fechner's law in
% metacognition: a quantitative model of working memory conifdence", by
% R van den Berg, AH Yoo, WJ Ma (Psych Rev, 2017).
%
% For questions, bug reports, etc, please email nronaldvdberg@gmail.com

function plot_data = plot_single_fit(subjidx,mflags)

% load data and best fitting parameters
fprintf('Subject S%d\n',subjidx);
sdata = read_data(subjidx);
fitpars = fit_model(subjidx,mflags);

% ========== compute summary stats for subject and model fit ========== 
N_vec = [3 6];
conf_vec = 0:5;
conf_vec_emp = 0:5;
conf_vec_fit = 0:5;
ntrials_per_N_fit = 1e6;
nbin_emp = 15;
nbin_fit = 51;
bin_emp = linspace(-pi,pi,nbin_emp+1);
bin_emp = bin_emp(1:end-1)+diff(bin_emp(1:2))/2;
bin_fit = linspace(-pi,pi,nbin_fit+1);
bin_fit = bin_fit(1:end-1)+diff(bin_fit(1:2))/2;
for ii=1:numel(N_vec)
    ntrials_per_N_emp = sum(sdata.N==N_vec(ii));
    fdata = gen_fake_data(mflags,fitpars,N_vec(ii),ntrials_per_N_fit);
    for jj=1:numel(conf_vec)
        % compute counts
        cnt_emp = hist(sdata.error_vec(sdata.N==N_vec(ii) & sdata.conf_vec==conf_vec(jj)),bin_emp);
        cnt_fit = hist(fdata.error_vec(fdata.conf_vec==conf_vec(jj)),bin_fit);
        cnt_fit = cnt_fit * ntrials_per_N_emp/ntrials_per_N_fit * diff(bin_emp(1:2))/diff(bin_fit(1:2));
                
        % store for group plotting
        plot_data.bin_emp = bin_emp;
        plot_data.bin_fit = bin_fit;
        plot_data.cnt_emp(ii,jj,:) = cnt_emp;
        plot_data.cnt_fit(ii,jj,:) = cnt_fit;
        
    end
    plot_data.avg_conf_emp(ii) = mean(sdata.conf_vec(sdata.N==N_vec(ii)));
    plot_data.avg_conf_fit(ii) = mean(fdata.conf_vec);    
end
for ii=1:numel(N_vec)
    % compute summary stats
    fdata = gen_fake_data(mflags,fitpars,N_vec(ii),1e5);
    binidx = interp1(conf_vec_fit,1:numel(conf_vec_fit),fdata.gamma_vec,'nearest','extrap');
    for jj=1:numel(conf_vec_emp)
        cvar_emp_cnt(jj) = sum(sdata.N==N_vec(ii) & sdata.conf_vec==conf_vec_emp(jj));
        [cvar_emp_mean(jj), cvar_emp_CI(jj,:)] = circ_var_bs(sdata.error_vec(sdata.N==N_vec(ii) & sdata.conf_vec==conf_vec_emp(jj))');
    end
    for jj=1:numel(conf_vec_fit)
        cvar_fit_cnt(jj) = sum(binidx==jj);
        [cvar_fit_mean(jj), cvar_fit_CI(jj,:)] = circ_var_bs(fdata.error_vec(binidx==jj)',cvar_emp_cnt(jj));
    end
    for jj=1:numel(conf_vec_emp)
        pconf_emp(jj) = sum(sdata.N==N_vec(ii) & sdata.conf_vec==conf_vec_emp(jj))/sum(sdata.N==N_vec(ii));
        pconf_fit(jj) = mean(fdata.conf_vec==conf_vec_emp(jj));
    end
    plot_data.conf_vec_emp = conf_vec_emp;
    plot_data.conf_vec_fit = conf_vec_fit;
    plot_data.cvar_emp(ii,:) = cvar_emp_mean;
    plot_data.cvar_fit(ii,:) = cvar_fit_mean;
    plot_data.cvar_emp_cnt(ii,:) = cvar_emp_cnt;
    plot_data.cvar_fit_cnt(ii,:) = cvar_fit_cnt;
    plot_data.pconf_emp(ii,:) = pconf_emp;
    plot_data.pconf_fit(ii,:) = pconf_fit;    
end

% ========== rest is plotting ========== 
if nargout>0
    return
end

% PANEL 1: error histogram as function of confidence
figure
set(gcf,'Position',get(gcf,'Position').*[0.1 0.7 2 .75]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.7 2 .75]);
N_vec = [3 6];
conf_vec = 0:5;
ntrials_per_N_fit = 1e6;
nbin_emp = 15;
nbin_fit = 51;
bin_emp = linspace(-pi,pi,nbin_emp+1);
bin_emp = bin_emp(1:end-1)+diff(bin_emp(1:2))/2;
bin_fit = linspace(-pi,pi,nbin_fit+1);
bin_fit = bin_fit(1:end-1)+diff(bin_fit(1:2))/2;
maxy = -Inf;
for ii=1:numel(N_vec)
    ntrials_per_N_emp = sum(sdata.N==N_vec(ii));
    fdata = gen_fake_data(mflags,fitpars,N_vec(ii),ntrials_per_N_fit);
    for jj=1:numel(conf_vec)
        % compute counts
        cnt_emp = hist(sdata.error_vec(sdata.N==N_vec(ii) & sdata.conf_vec==conf_vec(jj)),bin_emp);
        cnt_fit = hist(fdata.error_vec(fdata.conf_vec==conf_vec(jj)),bin_fit);
        cnt_fit = cnt_fit * ntrials_per_N_emp/ntrials_per_N_fit * diff(bin_emp(1:2))/diff(bin_fit(1:2));
        
        % plot
        subplot(2,6,(ii-1)*6+jj);
        plot(bin_emp,cnt_emp,'ko','markersize',4,'markerfacecolor','k');
        hold on
        plot(bin_fit,cnt_fit,'k','linewidth',1);
        maxy=max(max(maxy,max(cnt_emp)),max(cnt_fit));
        drawnow
        
        % set axis limits, labels, etc
        xlim([-pi pi]);
        set(gca,'Xtick',[-pi -pi/2 0 pi/2 pi],'XtickLabel',{'-pi','','0','','pi'});
        set(gca,'Ytick',0:30:500);        
        ylabel('trial count');
        xlabel('error');
        box off
    end
    plot_data.avg_conf_emp(ii) = mean(sdata.conf_vec(sdata.N==N_vec(ii)));
    plot_data.avg_conf_fit(ii) = mean(fdata.conf_vec);    
end
for ii=1:12
    subplot(2,6,ii)
    ylim([0 maxy*1.1]);
end
subplot(2,6,1)
title(sprintf('Subject #%d',subjidx));

% PANEL 2: p(conf) and cvar vs conf
figure
set(gcf,'Position',get(gcf,'Position').*[0.1 0.1 2 .5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.1 2 .5]);
N_vec = [3 6];
conf_vec_emp = 0:5;
conf_vec_fit = 0:5;
histmax = -Inf;
for ii=1:numel(N_vec)
    % compute summary stats
    fdata = gen_fake_data(mflags,fitpars,N_vec(ii),1e5);
    binidx = interp1(conf_vec_fit,1:numel(conf_vec_fit),fdata.gamma_vec,'nearest','extrap');
    for jj=1:numel(conf_vec_emp)
        cvar_emp_cnt(jj) = sum(sdata.N==N_vec(ii) & sdata.conf_vec==conf_vec_emp(jj));
        [cvar_emp_mean(jj), cvar_emp_CI(jj,:)] = circ_var_bs(sdata.error_vec(sdata.N==N_vec(ii) & sdata.conf_vec==conf_vec_emp(jj))');
    end
    for jj=1:numel(conf_vec_fit)
        cvar_fit_cnt(jj) = sum(binidx==jj);
        [cvar_fit_mean(jj), cvar_fit_CI(jj,:)] = circ_var_bs(fdata.error_vec(binidx==jj)',cvar_emp_cnt(jj));
    end
    for jj=1:numel(conf_vec_emp)
        pconf_emp(jj) = sum(sdata.N==N_vec(ii) & sdata.conf_vec==conf_vec_emp(jj))/sum(sdata.N==N_vec(ii));
        pconf_fit(jj) = mean(fdata.conf_vec==conf_vec_emp(jj));
    end
    histmax=max(max(histmax,max(pconf_emp)),max(pconf_fit));
    
    % plot histogram of conf - discrete
    subplot(1,4,ii);
    hold on
    bar(0:2:11    ,pconf_emp,.25,'FaceColor',[0 0 0]);
    hold on;
    bar(0.5:2:11.5,pconf_fit,.25,'FaceColor',[.7 .7 .7]);
    xlim([-0.5 11]);
    xlabel('Reported confidence');
    ylabel('Frequency');
    set(gca,'XTick',.25:2:10.25,'XTickLabel',0:5);
    set(gca,'Ytick',0:.2:1);
    box off    
    
    % plot cvar vs conf
    subplot(1,4,ii+2);
    hold on
    errorbar(conf_vec_emp,cvar_emp_mean,cvar_emp_CI(:,1),cvar_emp_CI(:,2),'ko','markerfacecolor','k','linestyle','none');
    valididx_fit = find(conf_vec_fit<=5 & ~isnan(cvar_fit_mean));
    plot(conf_vec_fit(valididx_fit),cvar_fit_mean(valididx_fit),'k-','linewidth',2);
    xlim([-0.5 5.5]);
    ylim([0 1])
    set(gca,'Xtick',0:5,'Ytick',0:.2:1);
    xlabel('Confidence');
    ylabel('Circular variance');
    box off
end
subplot(1,4,1);
ylim([0 histmax*1.1]);
title(sprintf('Subject #%d',subjidx));
subplot(1,4,2);
ylim([0 histmax*1.1]);


%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Compute bootstrapped 95% confidence intervals on circular variance  %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
function [cvar_m, cvar_CI] = circ_var_bs(Y,n)
if numel(Y)==0
    cvar_m = NaN;
    cvar_CI = [NaN NaN];
    return
end
nsamples = 1000;
if ~exist('n','var')
    Q=randi(numel(Y),numel(Y),nsamples);
else
    Q=randi(numel(Y),n,nsamples);    
end
cvar_vec = circ_var(Y(Q));
cvar_m=mean(cvar_vec);
Z = sort(normrnd(cvar_vec,1e-5));
idx = round([.025 .975]*numel(Z));
idx = max(idx,1);
cvar_CI(1) = cvar_m-Z(idx(1));
cvar_CI(2) = Z(idx(2))-cvar_m;


%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Generate synthetic data using one of the models from the factorial  %
% model space                                                         %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
function fdata = gen_fake_data(mflags,mpars,N,n_trials)

% unpack parameters
J1    = mpars(1);
power = mpars(2);
tau   = mpars(3);
sigma = mpars(4);
a     = mpars(5);
b     = mpars(6);
Kpar  = mpars(7);
Q     = mpars(8);

% compute mapping between J and kappa
kmap = [linspace(0,10,250) linspace(10.001,1e4,250)];
Jmap = kmap.*besseli(1,kmap,1)./besseli(0,kmap,1);

% draw #remembered items for each trial
if mflags(1)==1
    K_vec = ones(1,n_trials)*N; % remember all items
elseif mflags(1)==2
    K_vec = ones(1,n_trials)*min(Kpar,N); % remember min(K,N) items, with K a fixed integer
elseif mflags(1)==3
    K_vec = min(poissrnd(Kpar,1,n_trials),N); % remember min(K,N) items, with K drawn from Poisson pdf for each trial
end

% draw precision for each trial
if mflags(2)==1
    % Q chunks, assigned as evenly as possible to remembered items
    Q_low = floor(Q./K_vec);   % NB: Q_low and Q_high depend on K and can thus differ per trial
    Q_high = Q_low+1;
    p_high = mod(Q,K_vec)./K_vec;
    enc_high = rand(1,n_trials)<p_high;
    J_vec = zeros(1,n_trials);
    J_vec(~enc_high) = Q_low(~enc_high)*J1;
    J_vec(enc_high)  = Q_high(enc_high)*J1;
elseif mflags(2)==2
    % Q chunks, assigned randomly to remembered items
    Q_vec = binornd(Q,1./K_vec);
    J_vec = Q_vec*J1;    
else
    % Resource is a continuous quantity: J = J1*min(K,N)^power
    J_vec = J1*min(K_vec,N).^power;
end

% add variability 
if mflags(3)==2
    J_vec = gamrnd(J_vec/tau,tau);
end

% add guesses 
p_encode = K_vec./N; % probability that target is among remembered items
enc_vec = rand(1,numel(p_encode))<p_encode;
J_vec(~enc_vec) = 1e-10;

% compute corresponding kappa values
kappa_vec = interp1(Jmap,kmap,J_vec,'nearest');

% draw estimation errors
fdata.error_vec = circ_vmrnd2(0,kappa_vec);
       
% compute confidence, add noise
fdata.gamma_vec = a*log(J_vec)+b;
if sigma>1e-10
    fdata.gamma_vec = normrnd(fdata.gamma_vec,sigma);
end

% find confidence bin for each trial
fdata.conf_vec=interp1(0:5,0:5,fdata.gamma_vec,'nearest','extrap');


%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Function to load data for specific subject  %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
function data = read_data(subjidx)
load RademakerEtAl2012Data;
data.error_vec = all_data{subjidx}.error_vec;
data.conf_vec  = all_data{subjidx}.conf_vec;
data.N         = all_data{subjidx}.N;
data.uN        = all_data{subjidx}.uN;
