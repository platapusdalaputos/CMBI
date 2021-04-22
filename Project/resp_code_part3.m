%% Code for 2nd order Polynomial 

iTest = iM+1; % testing image start
iLim = 1500-iM;

%% Model Images poly 2

% Deform the source image with the Model transformations

Z = 1;
for i = 1:iLim
    
        [modImg2(Z), modDef2(Z), modDis2(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_p2(:,i), cpg2_p2(:,i), SDMS, source_image, images(:,i+iM));

        Z = Z+1;
end
%% Masking

% Deform the source mask image with the model transformations

% Variable name = Deformed Source Map Images (mDSMI)

Z = 1;
for i = 1:iLim
    
        [mDSMI2(Z), mDSMI_def2(Z), mDSMI_dis2(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_p2(:,i), cpg2_p2(:,i), SDMS, source_mask, images(:,i+iM));

        Z = Z+1;
end

%% Map the mask onto the deformed modelled images


for i = 1:iLim
    for j = 1:160
        for k = 1:160
            if mDSMI2(:,i).img(j,k) == 0
                modImg2(:,i).img(j,k) = 0;
                modDef2(:,i).img(j,k) = 0;                
            end      
        end 
    end
end

%% Plot arbitrary figures to see if its working


figure(2);
subplot(1,3,1)
dispNiiSlice(modImg2(:,200),"z",1)
title('Model Deformed Masked MRI Image')

subplot(1,3,2)
dispNiiSlice(mDSMI2(:,200),"z",1)
title('Model Deformed Mask Image')

subplot(1,3,3)
dispNiiSlice(modDef2(:,200),"z",1)
title('Model Deformed Mask Deformation Field Image')



%% Deformation Field Error Again

% Create a new struct with the same values as rDSMI for the for loop
DEF_def22 = rDSMI_def;

% Calculate the deformation field error for every voxel/pixel among the...
% ...1400 target imgs.

for i = 1:iLim
    for j = 1:160
        for k = 1:160
            % Registration Deformation field minus Model Deformation field
         DEF_def22(i).img(j,k) = regIm_def(:,i).img(j,k) - modDef2(:,i).img(j,k);  
         
        end 
    end
end

%% Making the deformation with just the registration mask
DEF_def2_mask2 = DEF_def22;

for i = 1:iLim
    for j = 1:160
        for k = 1:160
            if rDSMI_def(:,i).img(j,k) == 0
                DEF_def2_mask2(i).img(j,k) = 0;
            end      
        end 
    end
end
%%
figure(3);
subplot(1,2,1)
dispNiiSlice(DEF_def2_mask(:,25),"z",1,[-2 2])
colorbar;
subplot(1,2,2)
dispNiiSlice(DEF_def2_mask2(:,25),"z",1,[-2 2])
colorbar;
%%
L2_norm2 = DEF_def2_mask2;


for i = 1:iLim
    for j = 1:160
        for k = 1:160
            % Registration Deformation field minus Model Deformation field
            L2_norm2(i).img(j,k) = sqrt((regIm_def(:,i).img(j,k) - modDef2(:,i).img(j,k))^2);
            
            if rDSMI_def(:,i).img(j,k) == 0
                L2_norm2(i).img(j,k) = 0;
            end      
        end 
    end
end
%% Calculating the L2 norm of each image again AA

L2_norm_b = DEF_def2_mask2;

Aa1 = zeros(160,160,10);
Aa2 = zeros(160,160,10);
L2Aa = zeros(160,160,10);
for i = 1:iLim

%     % Registration Deformation field minus Model Deformation field
%     L2_norm_a(i).img(:,:,1,1,1) = regIm_def(:,i).img(:,:,1,1,1) - modDef(:,i).img(:,:,1,1,1);
%     L2_norm_a(i).img(:,:,1,1,2) = sqrt((regIm_def(:,i).img(:,:,1,1,2) - modDef(:,i).img(:,:,1,1,2))^2);   
    
    Aa1(:,:,i) = regIm_def(:,i).img(:,:,1,1,1);
    Aa2(:,:,i) = modDef2(:,i).img(:,:,1,1,1);
    
end
%% Calculating the L2 norm of each image again BB



Bb1 = zeros(160,160,10);
Bb2 = zeros(160,160,10);
L2Bb = zeros(160,160,10);
for i = 1:iLim

    
    Bb1(:,:,i) = regIm_def(:,i).img(:,:,1,1,2);
    Bb2(:,:,i) = modDef2(:,i).img(:,:,1,1,2);
    
end
%%

for i = 1:iLim
    for j = 1:160
        for k = 1:160
            L2Aa(j,k,i) = sqrt((Aa1(j,k,i) - (Aa2(j,k,i)))^2);
            L2Bb(j,k,i) = sqrt((Bb1(j,k,i) - (Bb2(j,k,i)))^2);
            
            
            if isnan(L2Aa(j,k,i))
                 L2Aa(j,k,i)=0;
            end
            if isnan(L2Bb(j,k,i))
                 L2Bb(j,k,i)=0;
            end
                
        end     
    end 
end

%% Calculating the mean
l2_mean_t = zeros(10,1);

for i = 1:iLim
    l2_mean_tAa(i,:) = mean(L2Aa(:,:,i),'all');
    l2_mean_tBb(i,:) = mean(L2Bb(:,:,i),'all');
end 


l2_mean_total_poly2 = mean(L2Aa+L2Bb,'all');

%% std
G = load('linear100trainmean.mat');
G1 = load('linear300trainmean.mat');
G2 = load('linear600trainmean.mat');
G3 = load('linear1000trainmean.mat');

gval = G.l2_mean_total;
gval1 = G1.l2_mean_total;
gval2 = G2.l2_mean_total;
gval3 = G3.l2_mean_total;

y = [gval, gval1, gval2, gval3];
x = [100, 300, 600, 1000];

G1 = load('poly100trainmean.mat');
G11 = load('poly300trainmean.mat');
G21 = load('poly600trainmean.mat');
G31 = load('poly21000trainmean.mat');

gval1 = G1.l2_mean_total_poly2;
gval11 = G11.l2_mean_total_poly2;
gval21 = G21.l2_mean_total_poly2;
gval31 = G31.l2_mean_total_poly2;

y1 = [gval1, gval11, gval21, gval31];
x1 = [100, 300, 600, 1000];

G12 = load('poly3_100trainmean.mat');
G112 = load('poly3_300trainmean.mat');
G212 = load('poly3_600trainmean.mat');
G312 = load('poly3_1000trainmean.mat');

gval12 = G12.l2_mean_total_poly3;
gval112 = G112.l2_mean_total_poly3;
gval212 = G212.l2_mean_total_poly3;
gval312 = G312.l2_mean_total_poly3;

y12 = [gval12, gval112, gval212, gval312];
x12 = [100, 300, 600, 1000];
figure;
plot(x,y,'-s','LineWidth',3)
hold on
plot(x1,y1,'-s','LineWidth',3)
hold on
plot(x12,y12,'-s','LineWidth',3)
legend('Linear','Poly 2^{nd}','Poly 3^{rd}')
xlabel('Training set, (No of images)','FontSize', 16)
ylabel('Average L2 norm, (\mu)','FontSize', 16)
xlim([0 1100])
ylim([1 3])

% 
% 
% 
% %%
% 
% mean_Val2 = [];
% 
% for i =1:1400
%     
%     SimpleArray2 = DEF_def2_mask2(:,i).img;
%     NewSimpleArray2 = SimpleArray2(SimpleArray2 ~= 0);
%     mean_Val2(i) = nanmean(NewSimpleArray2,'all');
%     
% end
% 
% %%
% 
% t2=1:iLim;
% figure(4);
% p2 = polyfit(t2,mean_Val2,1);
% f2 = polyval(p2,t2);
% % subplot(1,3,1)
% plot(t2, mean_Val2,'.',t2,f2,'-');
% hold on 
% plot(t2,f2,'r-');
% 
% %%
% figure(4);
% subplot(1,2,1)
% dispNiiSlice(L2_norm(:,9),"z",1,[-2 2])
% colorbar;
% subplot(1,2,2)
% dispNiiSlice(L2_norm2(:,9),"z",1,[-2 2])
% colorbar;