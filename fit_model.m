% function fitpars = fit_model(subjidx,mflags)
% 
% Generic function to fit one of the models from the factorial model
% space (see Table 3 of paper).
%
% The flags of the main model used in the paper (P-C-VP) are [3 3 2].
%
% INPUT
%  subjidx : subject index (integer in  range 1-6)
%  mflags  : 1x3 vector with model flags:
%            mflags(1): 1=All, 2=Fixed, 3=Poisson number of items encoded
%            mflags(2): 1=Quantized-even, 2=Quantized-random, 3=Continuous
%            mflags(3): 1=EP, 2=VP
% OUTPUT 
%  fitpars : maximum likelihood estimates of model parameters
%
% EXAMPLE
%  >> fit_model(1,[3 3 2])  % Fit main model to data of subject 1
%  
% This file is part of the code published with the paper "Fechner's law in
% metacognition: a quantitative model of working memory conifdence", by
% R van den Berg, AH Yoo, WJ Ma (Psych Rev, 2017).
%
% For questions, bug reports, etc, please email nronaldvdberg@gmail.com

function fitpars = fit_model(subjidx,mflags)

% check if result file exist; if so, load and return
results_fname = sprintf('saved_results/results_%d_%d_%d_%d.mat',subjidx,mflags);
if exist(results_fname,'file')
    load(results_fname,'maxLLH','AIC','fitpars');
    return
end

% check if we can use Matlab's parallel computing toolbox; if so, make sure
% that a pool is open (not essential, as it could wait until first real 
% parfor call, but doing it here avoids that ETL computation will be off)
parcomp = exist('parfor','builtin');
if parcomp  
    parfor ii=1,
        fprintf('Parallel Computing Toolbox opened\n');
    end
end

% Initialize structure with global variables. To speed up things, parameter 
% settings for the optimization method were lowered a bit compared to what 
% was used for the paper (see comments behind the lines below for paper
% settings)
gvar.nPar = 8;               % total number of parameters in the factorial space (J1, power, tau, sigma, c1, c5, Kpar, Q)
gvar.nMCSamples   = 500;    % number of MC samples to draw when computing model predictions (paper: 2500)
error_range = linspace(0,pi,91); % discretization of the error space
gvar.error_range = error_range(1:end-1)+diff(error_range(1:2))/2;
gvar.kappa_max      = 10000;
gvar.kappa_map      = [linspace(0,10,250) linspace(10.001,gvar.kappa_max,250)];
gvar.J_map          = gvar.kappa_map.*besseli(1,gvar.kappa_map,1)./besseli(0,gvar.kappa_map,1);
gvar.popSizeStart   = 1000;    % size of population at start (paper: 1500)
gvar.popSizeMin     = 100;     % minimum population size (paper: 200)
gvar.popSizeRate    = .98;     % size of new population is popSizeRate times size of current population (unless it is already at minimum)  (paper: .98)
gvar.nKids          = 1;       % number of kids per member (paper: 1)
gvar.nGenerations   = 250;     % number of generations to simulate (paper: 500)
gvar.ls = zeros(1,gvar.nPar);  % log sigma of the log normal distribution used to mutate parameter values
gvar.dls            = -0.01;   % difference in ls from generation to generation (paper: -0.01)
gvar.samplers = [2 NaN NaN 2 1 1 NaN NaN]; % 0=no mutation, 1=normal, 2=log normal, 3=uniform on [X-2, X+2]
gvar.mflags = mflags;        % add model flags for convenient passing to functions

% create matrix with initial parameter samples 
gvar.popSize = gvar.popSizeStart;
X_mat = zeros(gvar.popSize,gvar.nPar);
if mflags(2)<=2
    X_mat(:,1) = rand(gvar.popSize,1)*2 + .5;  % J1 / J1bar in discrete models (J = J1*#chunks)
    X_mat(:,2) = NaN(gvar.popSize,1);          % power is irrelevant in quantized models
    gvar.samplers(2) = 0;
elseif mflags(2)==3
    X_mat(:,1) = rand(gvar.popSize,1)*20 + 5;  % J1 / J1bar in continuous models (J = J1*N^power)
    X_mat(:,2) = rand(gvar.popSize,1)-2;       % power
    gvar.samplers(2) = 1;
end
if mflags(3)==1
    X_mat(:,3) = zeros(gvar.popSize,1);        % tau=0 (EP)
    gvar.samplers(3) = 0;
elseif mflags(3)==2
    X_mat(:,3) = rand(gvar.popSize,1)*10 + 1;  % tau>0 (VP)
    gvar.samplers(3) = 2;
end
X_mat(:,4) = rand(gvar.popSize,1)*1;           % sigma
X_mat(:,5) = rand(gvar.popSize,1)*2;           % a
X_mat(:,6) = rand(gvar.popSize,1)*4;           % b
if mflags(1)==1
    X_mat(:,7) = ones(gvar.popSize,1)*Inf;     % no slot limit
    gvar.samplers(7) = 0;
elseif mflags(1)==2
    X_mat(:,7) = randi(10,gvar.popSize,1);     % fixed slot limit
    gvar.samplers(7) = 3;
elseif mflags(1)==3
    X_mat(:,7) = rand(gvar.popSize,1)*10;      % mean of poisson distributed slot limit
    gvar.samplers(7) = 2;
end
if mflags(2)<=2
    X_mat(:,8) = randi(20,gvar.popSize,1);     % resource comes in Q discrete quanta
    gvar.samplers(8) = 3;
elseif mflags(2)==3
    X_mat(:,8) = Inf(gvar.popSize,1);          % resource is continuous
    gvar.samplers(8) = 0;
end

% compute total number of log likelihood (LLH) evaluations that will be performed
y(1)=gvar.popSizeStart;
for ii=2:gvar.nGenerations
    y(ii) = max(round(y(ii-1)*gvar.popSizeRate),gvar.popSizeMin);
end
n_eval_total = sum(y);

% load subject data and look up confidence and estimation error bin indices
data = read_data(subjidx);
for ii=1:length(data.uN)
    trial_idx = find(data.N==data.uN(ii));
    data.error_idx{ii} = interp1(gvar.error_range,1:length(gvar.error_range),abs(data.error_vec(trial_idx)),'nearest','extrap');
    data.conf_idx{ii} = interp1(0:5,1:6,data.conf_vec(trial_idx),'nearest','extrap');
end

% compute fitness of each member in initial generation 
fprintf('Computing fitness of first generation members...');
tic
LLH = zeros(1,gvar.popSize);
if parcomp
    parfor jj=1:gvar.popSize
        LLH(jj) = compute_LLH(X_mat(jj,:), data, gvar);
    end
else
    for jj=1:gvar.popSize
        LLH(jj) = compute_LLH(X_mat(jj,:), data, gvar);
    end
end
eval_cnt=gvar.popSize;
fprintf('\n');

% Run genetic algorithm
for ii=1:gvar.nGenerations

    % reproduce!
    X_mat_new = zeros(gvar.popSize*gvar.nKids,gvar.nPar);
    start_idx = 1;
    for jj=1:gvar.popSize        
        for kk=1:gvar.nPar
            X_mat_new(start_idx:start_idx+gvar.nKids-1,kk) = reproduce(X_mat(jj,kk),gvar.ls(kk),gvar.samplers(kk),gvar.nKids);
        end
        start_idx = start_idx + gvar.nKids;
    end
    
    % compute fitness of each member in new generation
    LLH_new = zeros(1,gvar.popSize*gvar.nKids);
    if parcomp
        parfor jj=1:gvar.popSize*gvar.nKids
            LLH_new(jj) = compute_LLH(X_mat_new(jj,:), data, gvar);
        end
    else
        for jj=1:gvar.popSize*gvar.nKids
            LLH_new(jj) = compute_LLH(X_mat_new(jj,:), data, gvar);
        end
    end
    eval_cnt=eval_cnt+numel(LLH_new);
    
    % merge old and new 
    X_mat_all = [X_mat; X_mat_new];
    LLH_all = [LLH LLH_new];

    % update population size
    gvar.popSize = max(gvar.popSizeMin,round(gvar.popSizeRate*gvar.popSize));

    % select fittest members (keep population size constant)
    [LLH_all, I] = sort(LLH_all,'descend');    
    LLH = LLH_all(1:gvar.popSize);
    X_mat = X_mat_all(I(1:gvar.popSize),:);
    fprintf('%d: mean/max LLH = %2.2f/%2.2f (ETL=%2.1f min)\n',ii,mean(LLH),max(LLH),toc/eval_cnt*(n_eval_total-eval_cnt)/60);
    fprintf('    %2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\n',X_mat(1,1),X_mat(1,2),X_mat(1,3),X_mat(1,4),X_mat(1,5),X_mat(1,6),X_mat(1,7),X_mat(1,8));
    
    % recompute LLH of current winner, to avoid pollution with "lucky" winners (LLH computation is stochastic)
    LLH(1) = compute_LLH(X_mat(1,:), data, gvar);
        
    gvar.ls = gvar.ls + gvar.dls;
end
        
% recompute LLHs 
fprintf('Recomputing LLH values with higher precision...');
for ii=1:numel(LLH)
    if parcomp
        parfor jj=1:20
            LLH_tmp(jj) = compute_LLH(X_mat(ii,:), data, gvar);
        end
    else
        for jj=1:20
            LLH_tmp(jj) = compute_LLH(X_mat(ii,:), data, gvar);
        end
    end
    LLH(ii) = mean(LLH_tmp);
end

% find ML parameter estimates
[maxLLH, max_idx] = max(LLH);
fitpars = X_mat(max_idx,:);

% compute BIC and estimate marginal model likelihood
nFreePars = sum(gvar.samplers~=0);
AIC = -2*maxLLH + 2*nFreePars;

% save results
save(results_fname,'fitpars','AIC','maxLLH','gvar');

%-%-%-%-%-%-%-%-%-%-%-%
% Likelihood function %
%-%-%-%-%-%-%-%-%-%-%-%
function LLH = compute_LLH(pars, data, gvar)

% unpack parameters
J1    = pars(1);
power = pars(2);
tau   = pars(3);
sigma = pars(4);
a     = pars(5);
b     = pars(6);
Kpar  = pars(7);
Q     = pars(8);

% compute log likelihood of data
LLH = 0;
for ii=1:length(data.uN)
    N = data.uN(ii);
    
    % draw number of *remembered* items for each MC sample
    if gvar.mflags(1)==1
        K_vec = ones(1,gvar.nMCSamples)*N;  % set number of slots equal to number of items -> all items remembered
    elseif gvar.mflags(1)==2
        K_vec = ones(1,gvar.nMCSamples)*min(Kpar,N); % fixed number of slots -> min(K,N) items remembered
    elseif gvar.mflags(1)==3
        K_vec = min(poissrnd(Kpar,1,gvar.nMCSamples),N); % Poisson distributed #slots -> min(K,N) items remembered, where K drawn from Poisson pdf in each MC sample
    end
    
    % draw J for each MC sample   
    if gvar.mflags(2)==1
        % Resource comes in Q discrete chunks that are distributed as evenly as possible across remembered items
        Q_low = floor(Q./K_vec);   % NB: Q_low and Q_high depend on K and can thus differ per trial
        Q_high = Q_low+1;
        p_high = mod(Q,K_vec)./K_vec;
        enc_high = rand(1,gvar.nMCSamples)<p_high;
        J_vec = zeros(1,gvar.nMCSamples);
        J_vec(~enc_high) = Q_low(~enc_high)*J1;
        J_vec(enc_high)  = Q_high(enc_high)*J1;
    elseif gvar.mflags(2)==2
        % Resource comes in discrete chunks that are distributed randomly across remembered items
        Q_vec = binornd(Q,1./K_vec);  % draw #chunks assigned to target (we have Q chunks and encode K_vec(i) items on each trial i)
        J_vec = Q_vec*J1;
    else
        % Resource is a continuous quantity: J = J1*min(K,N)^power
        J_vec = J1*min(K_vec,N).^power;
    end
    
    % add noise to J if it's a VP model
    if gvar.mflags(3)==2 && tau>0.001
        J_vec = gamrnd(J_vec/tau,tau);
    end
    
    % compute probability that target is among remembered items and stochastically insert guesses
    p_encode = K_vec./N;
    enc_vec = rand(1,numel(p_encode))<p_encode;
    J_vec(~enc_vec) = 0;  % set J to 0 on trials at which target was not among remembered items
        
    % make sure J is in numerically stable range
    J_vec = max(min(J_vec,max(gvar.J_map)),1e-6); 

    % compute corresponding kappas
    kappa = interp1(gvar.J_map,gvar.kappa_map,J_vec);    

    % compute probability of estimation error under each sample of J
    p_shat_mat = bsxfun(@rdivide,exp(bsxfun(@times,kappa',cos(gvar.error_range))),2*pi*besseli0_fast(kappa',0));

    % compute probability mass in each confidence bin under each sample of J
    gamma_vec = a*log(max(J_vec,1e-10))+b;
    cvec = [0.5 1.5 2.5 3.5 4.5]; % confidence criteria, in confidence space
    p_conf_mat = zeros(numel(J_vec),7);
    p_conf_mat(:,7) = 1;
    for jj=1:5
        p_conf_mat(:,jj+1) = normcdf(cvec(jj),gamma_vec,sigma);
    end
    p_conf_mat = diff(p_conf_mat,[],2);
    p_conf_mat(:,1) = p_conf_mat(:,1)+1e-20;  % to make sure that we don't have [0 0 0 0 0] for very small J. When J very small, we want to have [1 0 0 0 0] after the normalization in next line
    p_conf_mat = bsxfun(@rdivide,p_conf_mat',sum(p_conf_mat'))';     % normalize to make sure that sum equals 1
      
    % compute error pdf for each confidence level (compute joint under each sample of J, then marginalize)
    p_error = zeros(6,size(p_shat_mat,2));
    for jj=1:6
        p_error(jj,:) = mean(bsxfun(@times,p_shat_mat,p_conf_mat(:,jj)));
    end
    
    % read out response probability for each trial
    for jj=1:length(data.error_idx{ii})        
        p_resp = p_error(data.conf_idx{ii}(jj),data.error_idx{ii}(jj));
        LLH = LLH + log(p_resp);
    end    
end
% check if no NaNs (if so, set log likelihood to -Inf)
if isnan(LLH) || isinf(LLH)
    LLH=-Inf;
end

%-%-%-%-%-%-%-%-%-%-%
% Mutation function %
%-%-%-%-%-%-%-%-%-%-%
function X_new = reproduce(X_old,ls,samplerId,nKids)
if samplerId==0
    % no mutation
    X_new = ones(1,nKids)*X_old;
elseif samplerId==1
    % normal distribution
    X_new = normrnd(X_old,exp(ls),1,nKids);
elseif samplerId==2
    % log normal distribution
    X_new = exp(normrnd(log(X_old),exp(ls),1,nKids));
elseif samplerId==3
    % Uniform distribution on [X_old-1, X_old+1]
    X_new = max(X_old + randi(5,1,nKids)-3,0);
end

%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Function to load data for specific subject  %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
function data = read_data(subjidx)
load RademakerEtAl2012Data;
data.error_vec = all_data{subjidx}.error_vec;
data.conf_vec  = all_data{subjidx}.conf_vec;
data.N         = all_data{subjidx}.N;
data.uN        = all_data{subjidx}.uN;


% function I0 = besseli0_fast(kappa,scaleflag)
% 
% Returns bessel function for zeroth order, real arguments. Faster than the
% built in matlab function. If scaleflag is set to 1, the output is scaled
% by exp(kappa)
%
% Code copy-pasted from
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/129418,
% which apparently was based on the algorithm given in Numerical Recipes.
function I0 = besseli0_fast(x,scaleflag)
if ~exist('scaleflag','var')
    scaleflag=0;
end
ax = abs(x);
I0 = zeros(size(x));
% ax<3.75
k = find(ax<3.75);
y=x(k)./3.75;
y=y.^2;
I0(k)=1.0+y.*(3.5156229+y.*(3.0899424+y.*(1.2067492+y.*(0.2659732+y.*(0.360768e-1+y.*0.45813e-2)))));
if scaleflag
    I0(k)=I0(k)./exp(ax(k));
end
% ax>=3.75
k = find(ax>=3.75);
y=3.75./ax(k);
if scaleflag
    I0(k)=1./sqrt(ax(k)).*(0.39894228+y.*(0.1328592e-1+y.*(0.225319e-2+y.*(-0.157565e-2+y.*(0.916281e-2+y.*(-0.2057706e-1+y.*(0.2635537e-1+y.*(-0.1647633e-1+y.*0.392377e-2))))))));
else
    I0(k)=exp(ax(k))./sqrt(ax(k)).*(0.39894228+y.*(0.1328592e-1+y.*(0.225319e-2+y.*(-0.157565e-2+y.*(0.916281e-2+y.*(-0.2057706e-1+y.*(0.2635537e-1+y.*(-0.1647633e-1+y.*0.392377e-2))))))));
end