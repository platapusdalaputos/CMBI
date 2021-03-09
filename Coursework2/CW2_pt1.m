%% CW2 

% Sample length
SL = 20;

% Initialise the stochastic component with std = 0.2 and mean = 0
stan_d = 0.2;

% Mean for the stochastic component
mu = 0;
mu_s1 = 1.5;
mu_s2 = 2.0;



%% Simulate the sample data 
% set the sample for both groups to 20 with mean 1.5 and 2.0 

% Initialising sample 1
s1 = stan_d.*randn(SL,1) + mu_s1;

% Initialising sample 2
s2 = stan_d.*randn(SL,1) + mu_s2;


Y = [s2;s1];

% Fix the random seed to student ID
rng(20206994);

% Calculation of the stochastic component for twenty samples
E = stan_d.*randn(SL,1) + mu;



%% Computing the t-statistic of the sample 1 and sample 2

[h,p,ci,stats] = ttest2(s1,s2);

% calculating the stats values for both samples x1 and x2
statsx1 = [mean(s1) std(s1) var(s1)];
statsx2 = [mean(s2) std(s2) var(s2)];
statsE = [mean(E) std(E) var(E)];



% Design matrix

X = [repmat([0 1],20,1); repmat([1 0],20,1)];

% Computing the Rank of the design matrix X and hence the coloumn space

dimCX = rank(X);

% compute the projection of y onto the column space
% 
XX = X'*X;
PX = (X*pinv(XX))*X' ;
check = PX^2;

% PX = dot(X,X)/dot(X,X);
tracc = trace(PX);

% finding the value of y_hat
y_hat = PX*Y;

% identity matrix with the same values
I = eye(40,40);

% calculating Rx
Rx = I - PX;

% calculating e_hat
e_hat = Rx*Y;

% or you can also use
% e_hat2 = Y-y_hat;


% calculating the projection of Y on error space


% dimension of C(X)^{perp}
Dim_error = trace(Rx);


% angle between e_hat and Y_hat
CosTheta = max(min(dot(e_hat,y_hat)/(norm(e_hat)*norm(y_hat)),1),-1);
ThetaInDegrees = real(acosd(CosTheta));

% c.vii) Determine the parameter \Beta
Beta = pinv(XX)*X'*Y;

% c.viii) Variance of stochastic component
variance = (e_hat'*e_hat)/(40-2);

% c.ix) 

%Estimated covariance matrix
S_beta = variance*pinv(XX);

% estimated standard deviation
std_S_beta = sqrt(S_beta);


%% x. Null hypothesis

% brute force method

J = X.*[Beta(2),Beta(2)];
BF = J(:,1)+J(:,2);

BFBF = BF'*BF;
PBF = (BF*pinv(BFBF))*BF' ;
YBF = PBF*Y;
e_BF = Y-YBF;

% F-statistic

% 38 degrees of freedom

v1_BF = trace(PX-PBF);
v2_BF = trace(I-PX);


ABF = ((e_BF'*e_BF)-(e_hat'*e_hat))/v1_BF;
BBF = (e_hat'*e_hat)/v2_BF;
FBF = ABF/BBF;

% BF plot 

% figure(1);
% subplot(1,2,1)
% imagesc(BF);
% axis equal
% axis off
% colorbar
% title('BF')
% 
% subplot(1,2,2)
% imagesc(PBF);
% axis equal
% axis off
% colorbar
% title('PBF')

% same result is present.


%%
% Null hypothesis 
% contrast vector if the response is not affected by the stimulus 2
Cc = [1 -1]';

% Null vector
N = null(Cc');

% reduced model
X0 = (X*N);

% Additional error in reduced model X0
X0X0 = X0'*X0;
PX0 = (X0*pinv(X0X0))*X0' ;
e_X0 = (I - PX0)*Y;

% F-statistic

% 38 degrees of freedom 
v1 = trace(PX-PX0);
v2 = trace(I-PX);

% My Method
A = ((e_X0'*e_X0)-(e_hat'*e_hat))/v1;
B = (e_hat'*e_hat)/v2;
F = A/B;


% Gary's method 
SSR_X  = e_hat'*e_hat;
SSR_X0 = e_X0'*e_X0;
eDim = (length(Y)-rank(PX));
eDim_X0 = (length(Y)-rank(PX0));
varhat = SSR_X/eDim;
varhat_X0 = SSR_X0/eDim_X0;
eDimDiff = eDim_X0 - eDim;
eIncrease = (SSR_X0-SSR_X)/eDimDiff;

Fstat = eIncrease/varhat;

%xii

% beta_X0 = pinv(X0)*Y;

covHat = varhat*pinv(X'*X);
t = (Cc'* Beta)/(sqrt(Cc'*S_beta*Cc));

% Plotting the figures for X, XO, PX, PX0

% figure(2);
% subplot(1,2,1)
% imagesc(X);
% axis equal
% axis off
% colorbar
% title('X')
% 
% subplot(1,2,2)
% imagesc(PX);
% axis equal
% axis off
% colorbar
% title('PX')
% 
% figure(2);
% subplot(1,2,1)
% imagesc(X0);
% axis equal
% axis off
% colorbar
% title('X0')
% 
% subplot(1,2,2)
% imagesc(PX0);
% axis equal
% axis off
% colorbar
% title('PX0')



%% 1. d) Compute the t-statistic with a different model

% d.i) determine design matrix and column space
Xd = [repmat([1 0 1],20,1); repmat([1 1 0],20,1)];

% determining the rank and hence column space of the design matrix
% In linear algebra, the rank of a matrix A is the dimension of the vector 
% space generated (or spanned) by its columns. 


dimCXd = rank(Xd);

% d. ii)
XdXd = Xd'*Xd;
PXd = Xd*pinv(XdXd)*Xd' ;


% figure(2);
% imagesc(PXd);
% colorbar

% d.iii)

% Contrast vector
Cd = [0 1 -1]';

% Null vector
Nd = null(Cd');

% reduced model
X0_d = (Xd*Nd);

% d.iv) t-statistic

% beta values
Beta_d = pinv(XdXd)*Xd'*Y;

% calculating Rx
Rxd = I - PXd;

% calculating e_hat
e_hat_d = Rxd*Y;

variance_d = (e_hat_d'*e_hat_d)/(40-2);

%Estimated covariance matrix
S_beta_d = variance*pinv(XdXd);

% Calculate t-statistic
td = (Cd'* Beta_d)/(sqrt(Cd'*S_beta_d*Cd));

% d.v)

%% 1.e) Compute the t statistic with another model

% determine Design matrix 
Xe = [repmat([1 0],20,1); repmat([1 1],20,1)];

% determine dimension
dimXe = rank(Xe);


%% 2.a) Paired T-test

% using ttest

[T_paired, pv] = ttest(s1,s2);

% different pvalues for ttest and ttest2


% b.i)

% Produce design matrix of 
SI = fliplr(eye(40));

% Append the old x0b
X2a = [repmat([1 0],39,1); repmat([1 1],1,1)];

designX = [X2a SI];

% plot design matrix

figure(3);
imagesc(designX);

designX_Rank = rank(designX);




%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                             PART 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1.a) Simulate a new dataset n1 = 6, n2 = 8

% Initialising sample 1
n1 = randn(6,1);

% Initialising sample 2
n2 = randn(8,1);


[h_n,p_n,ci_n,stats_n] = ttest2(n1,n2);


% b.i) 

D = [n1;n2];

% b.ii)

%Compute all the permutations of D

c1 = combnk(D(1:6),1);
c2 = combnk(D(7:14),1);
permuteD = combvec(c1',c2')';

% compute the corresponding test statistics for each "group"
% Produce ordered group

% [H,PP,ci,sTatsg1] = ttest2(

% for i=1:length(permuteD)
%     
%    [H,PP,ci,sTats] = ttest2(permuteD(i,1),permuteD(i,2));
%     
% end




[H,PP,ci,sTats1] = ttest2(permuteD(1:6,1),permuteD(1:6,2));
[H,PP,ci,sTats2] = ttest2(permuteD(7:13,1),permuteD(7:13,2));
[H,PP,ci,sTats3] = ttest2(permuteD(14:20,1),permuteD(14:20,2));
[H,PP,ci,sTats4] = ttest2(permuteD(21:27,1),permuteD(21:27,2));
[H,PP,ci,sTats5] = ttest2(permuteD(28:34,1),permuteD(28:34,2));
[H,PP,ci,sTats6] = ttest2(permuteD(35:41,1),permuteD(35:41,2));
[H,PP,ci,sTats7] = ttest2(permuteD(42:48,1),permuteD(42:48,2));

% J = extractfield(sTats,'tstat');
%%
    
ggg = struct2cell(sTats1);
fff = cell2mat(ggg);

ggg1 = struct2cell(sTats2);
fff1 = cell2mat(ggg1);

ggg2 = struct2cell(sTats3);
fff2 = cell2mat(ggg2);

ggg3 = struct2cell(sTats4);
fff3 = cell2mat(ggg3);

ggg4 = struct2cell(sTats5);
fff4 = cell2mat(ggg4);

ggg5 = struct2cell(sTats6);
fff5 = cell2mat(ggg5);

ggg6 = struct2cell(sTats7);
fff6 = cell2mat(ggg6);

FF = [fff fff1 fff2 fff3 fff4 fff5 fff6];

figure(3);
histogram(FF,6);




