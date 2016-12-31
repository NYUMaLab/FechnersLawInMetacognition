% function create_Fig4A
% 
% Reproduces Fig 4A of the paper "Fechner's law in metacognition: a 
% quantitative model of working memory conifdence" (Van den Berg et al.
% 2017, Psych Rev)
%
% For questions, bug reports, etc, please email nronaldvdberg@gmail.com

function create_Fig4A
figure
set(gcf,'Position',get(gcf,'Position').*[0.1 0.1 0.5 0.5]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.1 0.5 0.5]);
hold on

lambda = [-1  -.5 -.25 .0001 .25 .5 .75 1];
J = linspace(.1,10,100);
plot(J,log(J),'k','linewidth',2);
for ii=1:numel(lambda)
    plot(J,(J.^lambda(ii)-1)./lambda(ii),'r--','linewidth',1);
end
ylim([0.1 5]);
xlabel('Memory precision, J');
ylabel('Memory confidence, \gamma');