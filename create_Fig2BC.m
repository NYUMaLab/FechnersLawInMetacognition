% function create_Fig2BC
% 
% Reproduces Fig 2B-C of the paper "Fechner's law in metacognition: a 
% quantitative model of working memory conifdence" (Van den Berg et al.
% 2017, Psych Rev)
%
% For questions, bug reports, etc, please email nronaldvdberg@gmail.com

function create_Fig2BC

% mapping between kappa and J
kmap = [linspace(0,10,250) linspace(10.001,500,250)];
Jmap = kmap.*besseli(1,kmap,1)./besseli(0,kmap,1);

% parameters for simulation
a = [1     1 1.5];
b = [3.5 3.5 3.5];
sigma_mc = [0.5 2.5];

% SIMULATE
figure
set(gcf,'Position',get(gcf,'Position').*[0.1 0.1 1 .5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.1 1 .5]);
subplot(1,2,1)
hold on
ntrials=50000;
J = gamrnd(.1,10,1,ntrials);
kappa = interp1(Jmap,kmap,J);
error = circ_vmrnd2(0,kappa);
conf1 = a(1)*log(J)+b(1) + randn(1,ntrials)*sigma_mc(1);
conf2 = a(2)*log(J)+b(2) + randn(1,ntrials)*sigma_mc(2);

% CREATE PANEL 2B
nplot = 1000;  % number of data points to include in scatter plot
Jplot = J(1:nplot);
confplot = conf1(1:nplot);
errorplot = error(1:nplot);
Jbinborders = [0 .5 5 Inf]; % border that define "low", "medium", "high" precision in Fig 2
cols = [.8 0 0; 0 .8 0; 0 0 .8];
syms = 'o+o';
for ii=1:numel(Jbinborders)-1
    idx = find(Jplot>Jbinborders(ii) & Jplot<Jbinborders(ii+1));
    if ii==1
        plot(confplot(idx),errorplot(idx),sprintf('r%s',syms(ii)),'color',cols(ii,:),'markersize',4);
    else
        plot(confplot(idx),errorplot(idx),sprintf('r%s',syms(ii)),'color',cols(ii,:),'markerfacecolor',cols(ii,:),'markersize',4);
    end
end 
xlabel('Confidence rating');
ylabel('Estimation error');
xlim([0 7]);
ylim([-pi pi]);
set(gca,'Ytick',[-pi -pi/2 0 pi/2 pi],'Yticklabel',{'-90','-45','0','45','90'});

% CREATE PANEL 2C
subplot(1,2,2)
confbin = 0:5;
dconf = 1;
for ii=1:numel(confbin)
    idx = conf1>confbin(ii)-dconf/2 & conf1<confbin(ii)+dconf/2;
    cvar1(ii) = circ_var(error(idx)');
    idx = conf2>confbin(ii)-dconf/2 & conf2<confbin(ii)+dconf/2;
    cvar2(ii) = circ_var(error(idx)');
end
plot(0:5,cvar1,'ko-','markerfacecolor','k','markersize',3)
hold on
plot(0:5,cvar2,'ro-','markerfacecolor','r','markersize',3)
xlim([-0.5 5.5]);
set(gca,'Xtick',0:5)
ylim([0 1]);
xlabel('Confidence rating');
ylabel('Circular variance of error distribution');