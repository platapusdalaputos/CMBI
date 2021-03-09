%% Part 2

% Initialise the stochastic component with std = 0.2 and mean = 0
stan_d = 0.2;

% Mean for the stochastic component
mu_n1 = 1.5;
mu_n2 = 2.0;

% Initialising sample 1
n1 = stan_d.*randn(6,1) + mu_n1;

% Initialising sample 2
n2 = stan_d.*randn(8,1) + mu_n2;

[h_n,p_n,ci_n,stats_n] = ttest2(n1,n2);

%% b.i) 

D = [n1;n2];

% b.ii)

%Compute all the permutations of D

c1 = combnk(D(1:6),1);
c2 = combnk(D(7:14),1);
permuteD = combvec(c1',c2')';

%% b.iii) compute the corresponding test statistics for each "group"

L=1;
hist = [];

for i=1:8
         
    [h, p, conf, stats] = ttest2(permuteD(L:L+5,1),permuteD(L:L+5,2));
     G = cell2mat(struct2cell(stats));
     hist(i) = G(1,:);

     L=L+6;

end

%% find the percentage of tstats greater than the original value

orig = cell2mat(struct2cell(stats_n));
tstat_orig = orig(1);

k = 0;
thresh =[];
for i=1:length(hist)
    J = hist(i);
    if J >= tstat_orig
        k=k+1;
        thresh(i) = J;
    end         
    
    perc = k/length(hist)*100;
end
%%
% b.iv) plot histogram with original tstat value

figure(1);

histogram(hist);
% hold on
% histogram(thresh);
hold on
xline(tstat_orig,'--r');
xlabel('tstat');
title('tstat of permutated D values');
legend('t-stat', 'p-value');
%% c) Mean differences instead of tstat


L1=1;
hist2 = [];

for i=1:8
         
    hist2(i) = mean(permuteD(L1:L1+5,1)) - mean(permuteD(L1:L1+5,2));

     L1=L1+6;

end

figure(2);

histogram(hist2);
% hold on
% histogram(thresh);
% hold on
% xline(tstat_orig,'--r');
xlabel('\mu_1 - \mu_2');
title('difference between means');
legend('\mu_1 - \mu_2', 'p-value');

%% d.i) 

noOfperms = 1000;
Dpermuted = zeros(length(D),noOfperms);
histd =  zeros(1,noOfperms);
for i =1:noOfperms
    
%     permutediD = randperm(length(D));
    Dpermuted(:,i) = randperm(length(D));
    D1 = D(Dpermuted(1:6,i));
    D2 = D(Dpermuted(7:14,i));
    
    [h1d, p1d, conf_d, stats_d] = ttest2(D1, D2);
     Gd = cell2mat(struct2cell(stats_d));
     histd(i) = Gd(1,:);
%      histd(i) = mean(D) - mean(Dpermuted(:,i));
    
    
end
%%
figure(3);

histogram(histd,10);




