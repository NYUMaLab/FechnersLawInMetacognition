% function gpd = plot_group_fit(subjidx,mflags)
% 
% Produces a plot with the subject-averaged maximum likelihood fit of the 
% model specified by mflags. The flags of the main model used in the 
% paper (P-C-VP) are [3 3 2].
%
% INPUT
%  mflags  : 1x3 vector with model flags:
%            mflags(1): 1=All, 2=Fixed, 3=Poisson number of items encoded
%            mflags(2): 1=Quantized-even, 2=Quantized-random, 3=Continuous
%            mflags(3): 1=EP, 2=VP
%
% OUTPOUT
%  gpd     : group plot data -- used in create_Fig8
%
% EXAMPLE
%  >> plot_group_fit([3 3 2])  % Plot fit of main model (P-C-VP)
%  
% This file is part of the code published with the paper "Fechner's law in
% metacognition: a quantitative model of working memory conifdence", by
% R van den Berg, AH Yoo, WJ Ma (Psych Rev, 2017).
%
% For questions, bug reports, etc, please email nronaldvdberg@gmail.com

function gpd = plot_group_fit(mflags)

% Get individual-subject plots
nSubj = 6; % number of subjects
for ii=1:nSubj
    plot_data = plot_single_fit(ii,mflags);
    bin_emp = plot_data.bin_emp;
    bin_fit = plot_data.bin_fit;
    cnt_emp(ii,:,:,:) = plot_data.cnt_emp;
    cnt_fit(ii,:,:,:) = plot_data.cnt_fit;
    conf_vec_emp = plot_data.conf_vec_emp;
    conf_vec_fit = plot_data.conf_vec_fit;
    cvar_emp(ii,:,:) = plot_data.cvar_emp;
    cvar_fit(ii,:,:) = plot_data.cvar_fit;
    cvar_emp_cnt(ii,:,:) = plot_data.cvar_emp_cnt; % number of data points that went into the cvar  
    cvar_fit_cnt(ii,:,:) = plot_data.cvar_fit_cnt; 
    pconf_emp(ii,:,:) = plot_data.pconf_emp;
    pconf_fit(ii,:,:) = plot_data.pconf_fit;
end

% compute weighted mean and std errors (weight by how many data points went into each cvar)
for ii=1:2 % loop over set sizes
    [cvar_m_emp(ii,:), cvar_s_emp(ii,:)] = weighted_mean_and_std(squeeze(cvar_emp(:,ii,:)),squeeze(cvar_emp_cnt(:,ii,:)));
    [cvar_m_fit(ii,:), cvar_s_fit(ii,:)] = weighted_mean_and_std(squeeze(cvar_fit(:,ii,:)),squeeze(cvar_fit_cnt(:,ii,:)));
    cvar_s_emp(ii,:) = cvar_s_emp(ii,:)/sqrt(size(cvar_emp,1));
    cvar_s_fit(ii,:) = cvar_s_fit(ii,:)/sqrt(size(cvar_fit,1));
end
% compute histogram means and std
hist_m_emp = squeeze(mean(pconf_emp,1));
hist_m_fit = squeeze(mean(pconf_fit,1));
hist_s_emp = squeeze(std(pconf_emp,1)/sqrt(6));
hist_s_fit = squeeze(std(pconf_fit,1)/sqrt(6));

% prepare output (used in create_Fig8)
gpd.cvar_m_emp = cvar_m_emp;
gpd.cvar_s_emp = cvar_s_emp;
gpd.hist_m_emp = hist_m_emp;
gpd.hist_s_emp = hist_s_emp;
gpd.cvar_m_fit = cvar_m_fit;
gpd.cvar_s_fit = cvar_s_fit;
gpd.hist_m_fit = hist_m_fit;
gpd.hist_s_fit = hist_s_fit;
gpd.conf_vec_fit = conf_vec_fit;
gpd.conf_vec_emp = conf_vec_emp;

if nargout>0
    return
end

% ====== error histogram as function of confidence ======
figure;
set(gcf,'Position',get(gcf,'Position').*[0.1 0.1 2 1]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.1 2 1]);
N_vec = [3 6];
conf_vec = 0:5;
maxy = -Inf;
for ii=1:numel(N_vec)
    for jj=1:numel(conf_vec)
        % compute mean and std for data and fit
        cnt_emp_mean = mean(squeeze(cnt_emp(:,ii,jj,:)),1);
        cnt_fit_mean = mean(squeeze(cnt_fit(:,ii,jj,:)),1);
        cnt_emp_eb = std(squeeze(cnt_emp(:,ii,jj,:)),[],1)/sqrt(size(cnt_emp,1));
        cnt_fit_eb = std(squeeze(cnt_fit(:,ii,jj,:)),[],1)/sqrt(size(cnt_fit,1));
        
        % plot
        subplot(2,6,(ii-1)*6+jj);
        hold on
        XX = [bin_fit bin_fit(end:-1:1)];
        YY = [cnt_fit_mean+cnt_fit_eb cnt_fit_mean(end:-1:1)-cnt_fit_eb(end:-1:1)];
        patch(XX,YY,[.7 .7 .7],'linestyle','none');
        errorbar(bin_emp,cnt_emp_mean,cnt_emp_eb,'ko','linestyle','none','markersize',4,'markerfacecolor','k');
        plot(bin_fit,cnt_fit_mean,'k','linewidth',1);
        maxy=max(max(maxy,max(cnt_emp_mean+cnt_emp_eb)),max(cnt_fit_mean+cnt_fit_eb));
    end
end

for ii=1:numel(N_vec)
    for jj=1:numel(conf_vec)
        subplot(2,6,(ii-1)*6+jj);
        xlim([-pi pi]);
        ylim([0 maxy*1.1]);
        set(gca,'Xtick',[-pi -pi/2 0 pi/2 pi],'XtickLabel',{'-pi','','0','','pi'});
        set(gca,'Ytick',[0 50 100]);
        if jj==1
            ylabel('trial count');
        else
            set(gca,'YtickLabel',{'','',''});
        end
        if ii==2
            xlabel('error');
        else
            set(gca,'XtickLabel',{'','','','',''});
        end
        title(sprintf('N=%d, conf=%d',N_vec(ii),conf_vec(jj)));
    end
end

% ======= histograms and cvar vs conf =========
figure;
set(gcf,'Position',get(gcf,'Position').*[0.1 0.1 2 .6]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.1 2 .6]);
N_vec = [3 6];
conf_vec_emp = 0:5;
for ii=1:numel(N_vec)    
    subplot(1,4,ii);
    hold on
    bar(0:3:15,hist_m_emp(ii,:),.3,'facecolor','k');
    bar(1:3:16,hist_m_fit(ii,:),.3,'facecolor',[.6 .6 .6]);
    errorbar(0:3:15,hist_m_emp(ii,:),hist_s_emp(ii,:),'k-','linestyle','none');
    errorbar(1:3:16,hist_m_fit(ii,:),hist_s_fit(ii,:),'k-','linestyle','none');
    xlim([-1 17])
    set(gca,'xtick',.5:3:15.5,'Xticklabel',0:5);
    xlabel('Confidence');
    ylabel('Proportion of trials');
    ylim([0 .5]);
end
for ii=1:numel(N_vec)
    % plot result
    subplot(1,4,ii+2);
    hold on
    XX = [conf_vec_fit conf_vec_fit(end:-1:1)];
    YY = [cvar_m_fit(ii,:)+cvar_s_fit(ii,:) cvar_m_fit(ii,end:-1:1)-cvar_s_fit(ii,end:-1:1)];
    patch(XX,YY,[.7 .7 .7],'linestyle','none');
    errorbar(conf_vec_emp,cvar_m_emp(ii,:),cvar_s_emp(ii,:),'ko-','linestyle','none','markersize',8,'markerfacecolor','k');
    plot(conf_vec_fit,cvar_m_fit(ii,:),'k','linewidth',1);
    xlim([-0.5 5.5]);
    ylim([0 1])
    set(gca,'Xtick',0:5,'Ytick',0:.2:1);
    xlabel('Confidence');
    ylabel('Circular variance');
end


function [m, s] = weighted_mean_and_std(A,cnt)

if ~all(size(A)==size(cnt))
    error('A and cnt should have equal dimensions');
end

if isvector(A)
    w = cnt./sum(cnt);
    m = nansum(w(w>0).*A(w>0));
    s = sqrt(nansum(w(w>0).*(A(w>0)-m).^2));
elseif ismatrix(A)
    for ii=1:size(A,2)
        [m(ii), s(ii)]= weighted_mean_and_std(A(:,ii),cnt(:,ii));
    end
end