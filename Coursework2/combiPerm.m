%% Function for calculating ttest for membership permutations 

function [ttestss, meanT] = combiPerm(n1,n2)

% Find the combinations of groups with numerical value b1 and/or b2 
% To take in consideration n1 = b1 and n2 = b2
% Then using set diff to find the values which aren't present


D = [n1;n2];

% numerical values 1 to 14 to represent listing in D

Dlist = 1:length(D);

% All possible combination from 1:14 in groups of 8
perM = combnk(Dlist, length(n1));
perM2 = zeros(length(perM),length(n2));

ttestss = zeros(length(perM), 1);

% Calculating the values in Dlist which are not in perM

for i = 1:length(perM)
    
    perM2(i,:) = setdiff(Dlist, perM(i,:));
    
end

% Mapping the permutation onto the actual values of 

permD = D(perM);
permD2 = D(perM2);
meanT = zeros(length(perM), 1);


for i = 1:length(perM)
    
     [~, ~, ~, stats]= ttest2(permD(i,:), permD2(i,:));
  
     G = cell2mat(struct2cell(stats));
     
     ttestss(i) = G(1,:);
     
     meanT(i) = mean(permD(i,:)) - mean(permD2(i,:));
        
end 

end