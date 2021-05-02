function [cpg1_lin, cpg2_lin, cpg1_p2, cpg2_p2, cpg1_p3, cpg2_p3] = cpgRetrieve(iLim, iTest, cpg1, cpg2, x_value, C, C_p2, C_p3)

% Copy the header from all of region 1 and region 2

Z = 1;
for i = 1:iLim
    hdr1(Z) = cpg1(:,i).hdr;
    hdr2(Z) = cpg2(:,i).hdr;
    Z = Z+1;
end 

% Retrieve the signal for the 1400 images and input into linear model

% Re-label the variables and matrices to understand whats going on

% linear
linearsurrogateSignals = [x_value(iTest:1500),ones(length(x_value(iTest:1500)),1)];
linearcoefficients = C;
lineartransformations = linearsurrogateSignals*linearcoefficients;

% retrieve cpg's
[cpg1_lin, cpg2_lin] = transForm(lineartransformations,hdr1,hdr2,iLim);

% Polynomial of 2nd order

p2surrogateSignals = [x_value(iTest:1500).^2,x_value(iTest:1500),ones(length(x_value(iTest:1500)),1)];
p2coefficients = C_p2;
p2transformations = p2surrogateSignals*p2coefficients;

% retrieve cpg's
[cpg1_p2, cpg2_p2] = transForm(p2transformations,hdr1,hdr2,iLim);


% Polynomial of 3rd order

p3surrogateSignals = [x_value(iTest:1500).^3, x_value(iTest:1500).^2,x_value(iTest:1500),ones(length(x_value(iTest:1500)),1)];
p3coefficient = C_p3;
p3transformations = p3surrogateSignals*p3coefficient;

% retrieve cpg's
[cpg1_p3, cpg2_p3] = transForm(p3transformations,hdr1,hdr2,iLim);

