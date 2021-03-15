%% Part 2

%% 2.1a) Simulate data set

% Initialise the stochastic component with std = 0.2 and mean = 0
stan_d = 0.2;

% Mean for the stochastic component
mu_n1 = 1.5;
mu_n2 = 2.0;

% Initialising sample 1
n1 = stan_d.*randn(6,1) + mu_n1;

% Initialising sample 2
n2 = stan_d.*randn(8,1) + mu_n2;

[~,pval,~,stats] = ttest2(n2,n1);
G = cell2mat(struct2cell(stats));
tval = G(1,:);

%% 2.1b.i) - 2.1b.iv)

% Construct 1D array
D = [n2;n1];

% all valid permutations
[ttestt, meanT] = combiPerm(n2,n1);


% find the percentage of tstats greater than the original value

k = 0;
thresh =[];
for i=1:length(ttestt)
    J = ttestt(i);
    if J >= tval
        k=k+1;
        thresh(i) = J;
    end         
    
    perc = k/length(ttestt);
end

% plot the figure and display the p-values
figure(1);

histogram(ttestt,100,'FaceColor',[0.8500 0.3250 0.0980]);
xlabel('t-statistic');
legend('t-statistic');


%% c) Difference in means as test statistic

% estimating the pvalue

k = 0;

for i=1:length(meanT)
    J = meanT(i);
    if J >= tval
        k=k+1;
    end             
    perc_mT = k/length(meanT);
end

%
% Plot the difference in means

figure(2);

histogram(meanT,100,'FaceColor',[0 0.4470 0.7410]);
xlabel('t-statistic (\mu_{permD} - \mu_{permD2})');
legend('\mu_{permD} - \mu_{permD2}');

%% 2d.i) estimating p-value with tstat and 1000 permutation

noOfperms = 1000;
Dpermuted = zeros(length(D),noOfperms);
histd =  zeros(1,noOfperms);

for i =1:noOfperms
    
    Dpermuted(:,i) = randperm(length(D));
    D1 = D(Dpermuted(1:6,i));
    D2 = D(Dpermuted(7:14,i));
    
    [h1d, p1d, conf_d, stats_d] = ttest2(D1, D2);
     Gd = cell2mat(struct2cell(stats_d));
     histd(i) = Gd(1,:);
    
    
end

k = 0;
thresh =[];
for i=1:noOfperms
    J = histd(i);
    if J >= tval
        k=k+1;
    end         
    
    perc_iter = k/length(histd);
end

figure(3);

histogram(histd,100,'FaceColor',[0.9290 0.6940 0.1250]);
hold on
xlabel('t-statistic');
legend('1000-permutations');

%% 2d.iii Duplicates in permutations

% finding the number of duplicates within the array histd

counts=[];
 c = unique(histd); % the unique values in histd  
 for i = 1:length(c)
   counts(i,1) = sum(histd==c(i)); % number of times each unique value is repeated
 end
 
duplicates = c(counts>1);

%%  Q2.2a) import and convert the files 

g1 = ["CPA4_diffeo_fa.img" "CPA5_diffeo_fa.img" "CPA6_diffeo_fa.img" "CPA7_diffeo_fa.img" "CPA8_diffeo_fa.img" "CPA9_diffeo_fa.img" "CPA10_diffeo_fa.img" "CPA11_diffeo_fa.img"];
g2 = ["PPA3_diffeo_fa.img" "PPA6_diffeo_fa.img" "PPA9_diffeo_fa.img" "PPA10_diffeo_fa.img" "PPA13_diffeo_fa.img" "PPA14_diffeo_fa.img" "PPA15_diffeo_fa.img" "PPA16_diffeo_fa.img"];
wm = 'wm_mask.img';

[MaxStat, tstattt] = maxtstat(g1,g2,wm);

NewSimpleArray = tstattt(tstattt ~= 0);
%%

figure(4);
histogram(NewSimpleArray,100);
hold on
xline(MaxStat,'--r');
xlabel('t-statistic');
legend('t-statistic', 'Max Tstat');

%% 2.2b)

D = [g1';g2'];

% numerical values 1 to 14 to represent listing in D

Dlist = 1:length(D);

% All possible combination from 1:16 in groups of 8
perM = combnk(Dlist, length(g1));
perM2 = zeros(length(perM),length(g2));

ttestss = zeros(length(perM), 1);

% Calculating the values in Dlist which are not in perM

for i = 1:length(perM)
    
    perM2(i,:) = setdiff(Dlist, perM(i,:));
    
end

% Mapping the permutation onto the actual values of the 

permD = D(perM);
permD2 = D(perM2);


%% 
statistics = [];

for i = 1:length(permD)
    
    w1 = permD(i,:);
    w2 = permD2(i,:);
    
    [statistics(i),~] = maxtstat(w1,w2,wm);
%     [A(:,:,i), B(:,:,i)] = convertImage(w1,w2,wm);
    
    
end 
%%

figure;

histogram(statistics,100);

%% Dlist = 1:14;

% All possible combination from 1:14 in groups of 8
perM = combnk(Dlist, 8);

% size of the amount of combinations
Np = length(perM);

% 
perM2 = zeros(Np,6);
ttestss = zeros(Np, 1);
meanTtest = zeros(Np, 1);


perMM = zeros(Np,6);

for i = 1:Np
    
    perMM(i,:) = setdiff(Dlist, perM(i,:));
    
end

D1 = D(perM);
D2 = D(perMM);

for i = 1:Np
    
     [~, ~, ~, stats]= ttest2(D1(i,:), D2(i,:));
  
     G = cell2mat(struct2cell(stats));
     
     ttestss(i) = G(1,:);
     
%      mean(D(L1:L1+5,1)) - mean(permuteD(L1:L1+5,2));
        
end 