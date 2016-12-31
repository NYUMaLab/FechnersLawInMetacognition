% function create_Fig8
% 
% Reproduces Fig 8 of the paper "Fechner's law in metacognition: a 
% quantitative model of working memory conifdence" (Van den Berg et al.
% 2017, Psych Rev)
%
% For questions, bug reports, etc, please email nronaldvdberg@gmail.com

function create_Fig8

figure;
set(gcf,'Position',get(gcf,'Position').*[0.1 0.1 1.7 2]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.1 1.7 2]);

mflags = [1 3 1; 2 3 1; 2 2 1; 1 3 2];
mnames = {'All-Continuous-EqualPrecision (Palmer, 1990)','Fixed-Continuous-EqualPrecision (Zhang and Luck, 2008)','Fixed-QuantizedUnequal-EqualPrecision (Zhang and Luck, 2008)','All-Continuous-VariablePrecision (Van den Berg et al. 2012)'};

for ii=1:4
    plot_data = plot_group_fit(mflags(ii,:));    
    N_vec = [3 6];
    for jj=1:numel(N_vec)
        subplot(4,4,(ii-1)*4+jj);
        hold on
        bar(0:3:15,plot_data.hist_m_emp(jj,:),.3,'facecolor','k');
        bar(1:3:16,plot_data.hist_m_fit(jj,:),.3,'facecolor',[.6 .6 .6]);
        errorbar(0:3:15,plot_data.hist_m_emp(jj,:),plot_data.hist_s_emp(jj,:),'k-','linestyle','none');
        errorbar(1:3:16,plot_data.hist_m_fit(jj,:),plot_data.hist_s_fit(jj,:),'k-','linestyle','none');
        xlim([-1 17])
        set(gca,'xtick',.5:3:15.5,'Xticklabel',0:5);
        xlabel('Confidence');
        ylabel('Proportion of trials');
        ylim([0 .5]);
    end
    title(mnames{ii});     
    for jj=1:numel(N_vec)
        % plot result
        subplot(4,4,(ii-1)*4+jj+2);
        hold on
        XX = [plot_data.conf_vec_fit plot_data.conf_vec_fit(end:-1:1)];
        YY = [plot_data.cvar_m_fit(jj,:)+plot_data.cvar_s_fit(jj,:) plot_data.cvar_m_fit(jj,end:-1:1)-plot_data.cvar_s_fit(jj,end:-1:1)];
        patch(XX,YY,[.7 .7 .7],'linestyle','none');
        errorbar(plot_data.conf_vec_emp,plot_data.cvar_m_emp(jj,:),plot_data.cvar_s_emp(jj,:),'ko-','linestyle','none','markersize',3,'markerfacecolor','k');
        plot(plot_data.conf_vec_fit,plot_data.cvar_m_fit(jj,:),'k','linewidth',1);
        xlim([-0.5 5.5]);
        ylim([0 1])
        set(gca,'Xtick',0:5,'Ytick',0:.2:1);
        xlabel('Confidence');
        ylabel('Circular variance');
    end
    drawnow
end