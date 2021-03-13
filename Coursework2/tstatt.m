%% Function  for calculating the tstat

function t = tstatt(y1,y2)

Y = [y2;y1];
% Design matrix

n = length(Y);

X = [repmat([1 0 1],length(y1),1); repmat([1 1 0],length(y1),1)];

% compute the projection of y onto the column space
 
XX = X'*X;
PX = (X*pinv(XX))*X' ;

% identity matrix with the same values
I = eye(length(Y),length(Y));

% calculating Rx
Rx = I - PX;

% calculating e_hat
e_hat = Rx*Y;


% c.vii) Determine the parameter \Beta
Beta = pinv(XX)*X'*Y;

% c.viii) Variance of stochastic component
variance = (e_hat'*e_hat)/(n-rank(X));


%Estimated covariance matrix
S_beta = variance*pinv(XX);

% contrast vector if the response is not affected by the stimulus 2
Cc = [0 1 -1]';

% t statistic
t = (Cc'* Beta)/(sqrt(Cc'*S_beta*Cc));



end 