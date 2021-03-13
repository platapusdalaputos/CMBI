%% Part 2

%% 2.1

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

%% b.ii) 

D = [n1;n2];


%% B.iii 

% Find the combinations of groups of 6 or 8
% To take in consideration n1 = 8 or n2 = 6
% Then using set diff to find the values which aren't present


% numerical values 1 to 14 to represent listing in D

Dlist = 1:14;

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
        
end 

figure;
histogram(ttestss);


%% b.iii) compute the corresponding test statistics for each "group"

% L=1;
% hist = [];
% 
% for i=1:8
%          
%     [~, ~, ~, stats] = ttest2(permuteD(L:L+5,1),permuteD(L:L+5,2));
%      G = cell2mat(struct2cell(stats));
%      hist(i) = G(1,:);
% 
%      L=L+6;
% 
% end

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


%% Q2 

% 2.a)

% All the matching values

% img1b = 'CPA5_diffeo_fa.img';
% img1c = 'CPA6_diffeo_fa.img';
% img1d = 'CPA7_diffeo_fa.img';
% img1e = 'CPA8_diffeo_fa.img';
% img1f = 'CPA9_diffeo_fa.img';
% img1g = 'CPA10_diffeo_fa.img';
% img1h = 'CPA11_diffeo_fa.img';

% img2b = 'PPA6_diffeo_fa.img';
% img2c = 'PPA9_diffeo_fa.img';
% img2d = 'PPA10_diffeo_fa.img';
% img2e = 'PPA13_diffeo_fa.img';
% img2f = 'PPA14_diffeo_fa.img';
% img2g = 'PPA15_diffeo_fa.img';
% img2h = 'PPA16_diffeo_fa.img';

% g1 = ["CPA4_diffeo_fa.img" "CPA5_diffeo_fa.img" "CPA6_diffeo_fa.img" "CPA7_diffeo_fa.img" "CPA8_diffeo_fa.img" "CPA9_diffeo_fa.img" "CPA10_diffeo_fa.img" "CPA11_diffeo_fa.img"];
% g2 = ["PPA3_diffeo_fa.img" "PPA6_diffeo_fa.img" "PPA9_diffeo_fa.img" "PPA10_diffeo_fa.img" "PPA13_diffeo_fa.img" "PPA14_diffeo_fa.img" "PPA15_diffeo_fa.img" "PPA16_diffeo_fa.img"];
% g3 =  ["PPA14_diffeo_fa.img" "PPA15_diffeo_fa.img" "PPA16_diffeo_fa.img"];




% maxT1 = maxtstat(img1a,img1b,'wm_mask.img');
% maxT2 = maxtstat(img2a,img2b,'wm_mask.img');

%%  Q2 import and convert the files 

g1 = ["CPA4_diffeo_fa.img" "CPA5_diffeo_fa.img" "CPA6_diffeo_fa.img" "CPA7_diffeo_fa.img" "CPA8_diffeo_fa.img" "CPA9_diffeo_fa.img" "CPA10_diffeo_fa.img" "CPA11_diffeo_fa.img"];
g2 = ["PPA3_diffeo_fa.img" "PPA6_diffeo_fa.img" "PPA9_diffeo_fa.img" "PPA10_diffeo_fa.img" "PPA13_diffeo_fa.img" "PPA14_diffeo_fa.img" "PPA15_diffeo_fa.img" "PPA16_diffeo_fa.img"];


% data = zeros(8,40,40,40);
% data1 = zeros(8,40,40,40);

data = zeros(40,40,40,8);
data1 = zeros(40,40,40,8);

N = length(g1);

for i = 1:N
    
    fid = fopen(g1(i), 'r', 'l'); % little-endian
    H = fread(fid, 'float'); % 16-bit floating point
    data(:,:,:,i) = reshape(H, [40 40 40]); % dimension 40x40x40
    
    fid1 = fopen(g2(i), 'r', 'l'); % little-endian
    H1 = fread(fid1, 'float'); % 16-bit floating point
    data1(:,:,:,i) = reshape(H1, [40 40 40]); % dimension 40x40x40
    
end

fid2 = fopen('wm_mask.img', 'r', 'l'); % little-endian
data2 = fread(fid2, 'float'); % 16-bit floating point
wm_mask = reshape(data2, [40 40 40]); % dimension 40x40x40

%%


% fid = fopen('CPA4_diffeo_fa.img', 'r', 'l'); % little-endian
% data = fread(fid, 'float'); % 16-bit floating point
% data = reshape(data, [40 40 40]); % dimension 40x40x40
% 
% 
% fid1 = fopen('wm_mask.img', 'r', 'l'); % little-endian
% data1 = fread(fid1, 'float'); % 16-bit floating point
% wm_mask = reshape(data1, [40 40 40]); % dimension 40x40x40
% 
% fid2 = fopen('PPA3_diffeo_fa.img', 'r', 'l'); % little-endian
% data2 = fread(fid2, 'float'); % 16-bit floating point
% data2 = reshape(data2, [40 40 40]); % dimension 40x40x40





% GG = zeros(40,40,40);
% CPA = data;
% PPA = data2;
%% CPA and PPA mapping

% Mapping the wm_mask onto the CPA and PPA images

CPA = data;
PPA = data1;

for l = 1:8
    for i=1:40
        for j=1:40
            for k=1:40 
                 if wm_mask(k,j,i) == 0

                     CPA(k,j,i,l) = 0;
                     PPA(k,j,i,l) = 0;
                     

                end
            end
        end    
    end 
end



%% PPA mapping

% Mapping the wm_mask onto the PPA images
% PPA mapping
% 
% for i=1:40
%     for j=1:40
%         for k=1:40 
%              if wm_mask(k,j,i) == 0
%                  PPA(k,j,i) = 0;
%                                 
%             end
%         end
%     end    
% end


%% Computing the tstat on every voxel

tstat_value = [];
ttest_tstat = [];
y1hodl = [];
y2hodl = [];


for i=1:40
    for j=1:40
        for k=1:40 
            for l =1:8

%                  y1hodl = CPA(k,j,i);
%                  y2hodl = PPA(k,j,i);
                    
%                  y1hodl(:,:,:,l) = CPA(k,j,i,l)
                 y1hodl(l) = CPA(k,j,i,l);
                 y2hodl(l) = PPA(k,j,i,l);

%                  y2hodl(:,:,:,l) = PPA(k,j,i,l)
                 
%                  Y = [y1hodl';y2hodl'];
                 

%              if k>8
                 transY1 = y1hodl';
                 transY2 = y2hodl';
%                  Y = [transY1;transY2]
%                  tstat_value = tstatt(y1hodl,y2hodl);
              [~,~,~,Sts] = ttest2(transY1,transY2);
                 tstat_value(k,j,i) = tstatt(transY1,transY2);
%                  if isnan(tstat_value(k,j,i))
%                      tstat_value(k,j,i)=0;
%              end
%              end
                
            end
        end
    end    
end 
%%

% NewSimpleArray = tstat_value(tstat_value ~= 0);
% MaxTstat = max(NewSimpleArray);
% 
% % Z = find(~tstat_value);
% % [x,~] = ind2sub(size(tstat_value),Z);
% % tstat_value(x,:) = [];
% 
% %%
% figure;
% % if tstat
% histogram(NewSimpleArray);
% 
% 
% % tstattt = tstatt(transY1, transY2);
% 
% % Y = [transY1,transY2];
% % % Design matrix
% % 
% % X = [repmat([0 0 1],40,1); repmat([0 1 0],40,1)];
% % 
% % % compute the projection of y onto the column space
% %  
% % XX = X'*X;
% % PX = (X*pinv(XX))*X' ;
