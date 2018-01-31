% function create_Fig9A
% 
% Reproduces Fig 9A of the paper "Fechner's law in metacognition: a 
% quantitative model of working memory conifdence" (Van den Berg et al.
% 2017, Psych Rev)
%
% For questions, bug reports, etc, please email nronaldvdberg@gmail.com
 
function create_Fig9A

% parameters
mu_target    = 2;
sigma_target = 1.25;
a         = 2;
b         = 3;
sigma_mc  = 1;
nConf = 20; % number of confidence values

% ======= plot 1: internal representation ========
x = linspace(-10,10,100);
ptarget_x = normpdf(x,mu_target,sigma_target); % gaussian distributions
pfoil_x = normpdf(x,0,1);
ptarget_x = ptarget_x/sum(ptarget_x(:))/diff(x(1:2)); % normalize
pfoil_x = pfoil_x/sum(pfoil_x(:))/diff(x(1:2));

% ======= plot 2: decision variable ========
dx = -log(sigma_target) - 1/2.*( ((x-mu_target).^2)./sigma_target.^2 - x.^2);

% ======= plot 3: distribution of d ========
nSamples = 1e7;
samples_foil = randn(1,nSamples);
samples_target = sigma_target.*randn(1,nSamples) + mu_target;
d_foil = -log(sigma_target) - 1/2.*( ((samples_foil-mu_target).^2)./sigma_target.^2 - samples_foil.^2);
d_target = -log(sigma_target) - 1/2.*( ((samples_target-mu_target).^2)./sigma_target.^2 - samples_target.^2);
centers_target = linspace(-10,10,100);
centers_foil = centers_target;
counts_target = hist(d_target,centers_target);
counts_foil = hist(d_foil,centers_foil);
counts_target = counts_target/sum(counts_target(:))/diff(centers_target(1:2));
counts_foil = counts_foil/sum(counts_foil(:))/diff(centers_foil(1:2));

% ======= plot 4: confidence ratings ========
% get indices of responses
idx_CR = d_foil < 0; % correct reject. foil, respond new
idx_FA = d_foil >=0; % false alarm. foil, respond old
idx_Miss = d_target < 0; % miss. target, respond new
idx_Hit = d_target >= 0; % hit. target. respond old

% confidence values
conf_new = a.*log(abs(d_foil)) + b; 
conf_old = a.*log(abs(d_target)) + b;
conf_new = normrnd(conf_new,sigma_mc); % adding metacognitive noise
conf_old = normrnd(conf_old,sigma_mc);
conf_new = round(conf_new); % rounding
conf_old = round(conf_old);
conf_new(conf_new < 1) = 1; % satisfying confidence boundaries
conf_old(conf_old < 1) = 1;
conf_new(conf_new > 10) = 10;
conf_old(conf_old > 10) = 10;

pold = nan(1,nConf); % probability respond "old" for each confidence value
pnew = nan(1,nConf); % probability respond "new" for each confidence value
for iconf = 1:10
    idx_confCR = idx_CR & (conf_new == iconf);
    idx_confFA = idx_FA & (conf_new == iconf); 
    idx_confMiss = idx_Miss & (conf_old == iconf);
    idx_confHit = idx_Hit & (conf_old == iconf);
    
    pnew(iconf+nConf/2) = sum(idx_confFA); 
    pnew(nConf/2 + 1 - iconf) = sum(idx_confCR);
    pold(iconf+nConf/2) = sum(idx_confHit);
    pold(nConf/2 + 1 - iconf) = sum(idx_confMiss);
end
pnew = pnew./sum(pnew);
pold = pold./sum(pold);

% calculate confidence boundaries
confBounds = exp(([1.5:9.5]-b)./a);
confBounds = [-confBounds 0 confBounds];

% ======== plot ===========
color_foil = [.4 .5 .8]; % foil color
color_target = [.7 .6 .4]; % target color

figure
subplot(2,2,1); hold on
plot(x,ptarget_x,'Color',color_target); 
plot(x,pfoil_x,'Color',color_foil)
title('Observation')
xlim([x(1) x(end)])
ylabel('Probability')
ylabel('x')
hold off

subplot(2,2,2); hold on
plot(x,dx,'k')
plot([x(1) x(end)],[0 0],'--','Color',0.7*ones(1,3))
xlim([x(1) x(end)])
ylabel('d')
xlabel('x')
title('Decision variable')

subplot(2,2,3); hold on
plot([confBounds; confBounds],[zeros(1,length(confBounds)); 0.4*ones(1,length(confBounds))],'Color',0.7*ones(1,3));
plot(centers_foil,counts_foil,'Color',color_foil); hold on
plot(centers_target,counts_target,'Color',color_target);
xlim([-7 7])
ylabel('probability')
xlabel('d')
title('Distribution of d')

subplot(2,2,4); hold on
plot(1:nConf,pnew,'-','Color',color_foil);
plot(1:nConf,pold,'-','Color',color_target);
set(gca,'Xtick',[1:10 11:20],'XtickLabel',{10,[],[],[],6,[],[],[],2,[],[],2,[],[],[],6,[],[],[],10});
xlim([0 21])
xlabel('Confidence')
title('confidence distributions')
hold off