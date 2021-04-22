%% 6.5 Create the Masks and Transformations and Deformation Field Errors

% Pre-requisites

% Signed Distance Map Segmentation = SDMS
SDMS = load_untouch_nii('0007_sdt.nii.gz');

% Source mask
source_mask = load_untouch_nii('0007_mask.nii.gz');

% Source image
source_image = load_untouch_nii('0007.nii.gz');

iTest = iM+1; % testing image start
iLim = 1500-iM;
%% Registration Images

% Deform the source image with the registration transformations
% Registation images = regIm
Z = 1;
for i = iTest:1500
    
        [regIm(Z), regIm_def(Z), regIm_dis(Z)] = ...
            deformNiiWithCPGsSliding(cpg1(:,i), cpg2(:,i), SDMS, source_image, images(:,i));

        Z = Z+1;
end

%% Masking The Registration Images

% Deform the source mask image with the registration transformations
% The deformNiiWithCPGsSliding will also calculate the Deformation Field
% Signed Distance Map Segmentation = SDMS
% Source_image = SI
% Target_image = TI
% [The_Image , deformation_Field, distortion_Field] = ...
% ...deformNiiWithCPGsSliding(registration_Region1, registration_Region2, SDMS, SI, TI);

% Creating the registration Deformed Source Map Images (rDSMI)

Z = 1;
for i = iTest:1500
    
        [rDSMI(Z), rDSMI_def(Z), rDSMI_dis(Z)] = ...
            deformNiiWithCPGsSliding(cpg1(:,i), cpg2(:,i), SDMS, source_mask, images(:,i));

        Z = Z+1;
end



% Map the mask onto the registration images

%% Create New registration Masked images (NRMI)

NRMI = regIm;

for i = 1:iLim
    for j = 1:160
        for k = 1:160
            if rDSMI(:,i).img(j,k) == 0
                NRMI(:,i).img(j,k) = 0;
                rDSMI_def(:,i).img(j,k) = 0; 
                regIm_def(:,i).img(j,k) = 0; 
            end      
        end 
    end
end


%% Plot arbitrary figures to see if its working

% 
% figure(1);
% subplot(1,3,1)
% dispNiiSlice(NRMI(:,108),"z",1)
% title('Registration Deformed Masked MRI Image')
% 
% subplot(1,3,2)
% dispNiiSlice(regIm(:,108),"z",1)
% title('Registration MRI Image')
% 
% subplot(1,3,3)
% dispNiiSlice(rDSMI_def(:,108),"z",1)
% title('Registration Deformed Mask Deformation Field Image')

%% Copy the header from all of region 1 and region 2

Z = 1;
for i = 1:iLim
    hdr1(Z) = cpg1(:,i).hdr;
    hdr2(Z) = cpg2(:,i).hdr;
    Z = Z+1;
end 

%% Retrieve the signal for the 1400 images and input into linear model

% Re-label the variables and matrices to understand whats going on

% linear
linearsurrogateSignals = [x_20(iTest:1500),ones(length(x_20(iTest:1500)),1)];
linearcoefficients = C;
lineartransformations = linearsurrogateSignals*linearcoefficients;

%% retrieve cpg's
[cpg1_lin, cpg2_lin] = transForm(lineartransformations,hdr1,hdr2,iLim);

%% Polynomial of 2nd order

p2surrogateSignals = [x_20(iTest:1500).^2,x_20(iTest:1500),ones(length(x_20(iTest:1500)),1)];
p2coefficients = C_p2;
p2transformations = p2surrogateSignals*p2coefficients;

% retrieve cpg's
[cpg1_p2, cpg2_p2] = transForm(p2transformations,hdr1,hdr2,iLim);


%% Polynomial of 3rd order

p3surrogateSignals = [x_20(iTest:1500).^3, x_20(iTest:1500).^2,x_20(iTest:1500),ones(length(x_20(iTest:1500)),1)];
p3coefficient = C_p3;
p3transformations = p3surrogateSignals*p3coefficient;

% retrieve cpg's
[cpg1_p3, cpg2_p3] = transForm(p3transformations,hdr1,hdr2,iLim);
% %%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%                    % Linear Deformation Fielf Error %
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %% Model Images
% 
% % Deform the source image with the Model transformations
% 
% Z = 1;
% for i = 1:iLim
%     
%         [modImg(Z), modDef(Z), modDis(Z)] = ...
%             deformNiiWithCPGsSliding(cpg1_lin(:,i), cpg2_lin(:,i), SDMS, source_image, images(:,i+iM));
% 
%         Z = Z+1;
% end
% %% Masking
% 
% % Deform the source mask image with the model transformations
% 
% % Variable name = Deformed Source Map Images (mDSMI)
% 
% Z = 1;
% for i = 1:iLim
%     
%         [mDSMI(Z), mDSMI_def(Z), mDSMI_dis(Z)] = ...
%             deformNiiWithCPGsSliding(cpg1_lin(:,i), cpg2_lin(:,i), SDMS, source_mask, images(:,i+iM));
% 
%         Z = Z+1;
% end
% 
% %% Map the mask onto the deformed modelled images
% 
% 
% for i = 1:iLim
%     for j = 1:160
%         for k = 1:160
%             if mDSMI(:,i).img(j,k) == 0
%                 modImg(:,i).img(j,k) = 0;
%                 modDef(:,i).img(j,k) = 0;                
%             end      
%         end 
%     end
% end
% 
% %% Plot arbitrary figures to see if its working
% 
% 
% figure(2);
% subplot(1,3,1)
% dispNiiSlice(modImg(:,100),"z",1)
% title('Model Deformed Masked MRI Image')
% 
% subplot(1,3,2)
% dispNiiSlice(mDSMI(:,100),"z",1)
% title('Model Deformed Mask Image')
% 
% subplot(1,3,3)
% dispNiiSlice(modDef(:,100),"z",1)
% title('Model Deformed Mask Deformation Field Image')
% 
% %% Deformation Field Error
% 
% % Create a new struct with the same values as rDSMI for the for loop
% DEF_def = rDSMI_def;
% 
% % Calculate the deformation field error for every voxel/pixel among the...
% % ...1400 target imgs.
% 
% for i = 1:iLim
%     for j = 1:160
%         for k = 1:160
%             % Registration Deformation field minus Model Deformation field
%          DEF_def(i).img(j,k) = rDSMI_def(:,i).img(j,k) - modDef(:,i).img(j,k);  
%          
%         end 
%     end
% end
% 
% %% Deformation Field Error Again
% 
% % Create a new struct with the same values as rDSMI for the for loop
% DEF_def2 = rDSMI_def;
% 
% % Calculate the deformation field error for every voxel/pixel among the...
% % ...1400 target imgs.
% 
% for i = 1:iLim
%     for j = 1:160
%         for k = 1:160
%             % Registration Deformation field minus Model Deformation field
%          DEF_def2(i).img(j,k) = regIm_def(:,i).img(j,k) - modDef(:,i).img(j,k);  
%          
%         end 
%     end
% end
% %% Making the deformation with just the registration mask
% DEF_def2_mask = DEF_def2;
% 
% for i = 1:iLim
%     for j = 1:160
%         for k = 1:160
%             if rDSMI_def(:,i).img(j,k) == 0
%                 DEF_def2_mask(i).img(j,k) = 0;
%             end      
%         end 
%     end
% end
% 
% %%
%     
% figure;
% 
% dispNiiSlice(DEF_def2_mask(:,2),"z",1,[-2 2])
% colorbar;
% 
% 
% %% Calculating the L2 norm of each image
% 
% L2_norm = DEF_def;
% mean_Val = [];
% 
% 
% for i = 1:10
%     for j = 1:160
%         for k = 1:160
%             % Registration Deformation field minus Model Deformation field
%             L2_norm(i).img(j,k) = sqrt((regIm_def(:,i).img(j,k) - modDef(:,i).img(j,k))^2);
%             
%             if rDSMI_def(:,i).img(j,k) == 0
%                 L2_norm(i).img(j,k) = 0;
%             end      
%         end 
%     end
% end
% 
% %% Calculating the L2 norm of each image again AA
% 
% L2_norm_b = DEF_def;
% 
% Aa1 = zeros(160,160,10);
% Aa2 = zeros(160,160,10);
% L2Aa = zeros(160,160,10);
% for i = 1:iLim
% 
% %     % Registration Deformation field minus Model Deformation field
% %     L2_norm_a(i).img(:,:,1,1,1) = regIm_def(:,i).img(:,:,1,1,1) - modDef(:,i).img(:,:,1,1,1);
% %     L2_norm_a(i).img(:,:,1,1,2) = sqrt((regIm_def(:,i).img(:,:,1,1,2) - modDef(:,i).img(:,:,1,1,2))^2);   
%     
%     Aa1(:,:,i) = regIm_def(:,i).img(:,:,1,1,1);
%     Aa2(:,:,i) = modDef(:,i).img(:,:,1,1,1);
%     
% end
% %% Calculating the L2 norm of each image again BB
% 
% 
% 
% Bb1 = zeros(160,160,10);
% Bb2 = zeros(160,160,10);
% L2Bb = zeros(160,160,10);
% for i = 1:iLim
% 
%     
%     Bb1(:,:,i) = regIm_def(:,i).img(:,:,1,1,2);
%     Bb2(:,:,i) = modDef(:,i).img(:,:,1,1,2);
%     
% end
% %%
% 
% for i = 1:iLim
%     for j = 1:160
%         for k = 1:160
%             L2Aa(j,k,i) = sqrt((Aa1(j,k,i) - (Aa2(j,k,i)))^2);
%             L2Bb(j,k,i) = sqrt((Bb1(j,k,i) - (Bb2(j,k,i)))^2);
%             
%             
%             if isnan(L2Aa(j,k,i))
%                  L2Aa(j,k,i)=0;
%             end
%             if isnan(L2Bb(j,k,i))
%                  L2Bb(j,k,i)=0;
%             end
%                 
%         end     
%     end 
% end
% 
% %% Calculating the mean
% l2_mean_t = zeros(10,1);
% 
% for i = 1:iLim
%     l2_mean_tAa(i,:) = mean(L2Aa(:,:,i),'all');
%     l2_mean_tBb(i,:) = mean(L2Bb(:,:,i),'all');
% end 
% 
% 
% l2_mean_total = mean(L2Aa+L2Bb,'all');
% 
% %% std
% 
% G = load('linear100trainmean.mat');
% G1 = load('linear300trainmean.mat');
% G2 = load('linear600trainmean.mat');
% G3 = load('linear1000trainmean.mat');
% 
% gval = G.l2_mean_total;
% gval1 = G1.l2_mean_total;
% gval2 = G2.l2_mean_total;
% gval3 = G3.l2_mean_total;
% 
% y = [gval, gval1, gval2, gval3];
% x = [100, 300, 600, 1000];
% figure;
% plot(x,y,'-o')
% xlim([0 1500])
% ylim([0 3])
% %%
% figure;
% imshow(L2Aa(:,:,1));
% 
% 
% %%
% % %     Bb(i) = sqrt((regIm_def(:,i).img(:,:,1,1,2) - modDef(:,i).img(:,:,1,1,2))^2);
% %     
% %     for j = 1:160
% %         for k = 1:160
% % %              L2_norm_b(i).img(:,:,1,1,1) = sqrt(Aa1(j,k)^2 - Aa2(j,k)^2);
% %              HH(i) = Aa1(j,k) - Aa2(j,k);
% %             
% %             
% %         end
% %     end
% % 
% %     
% % 
% % end
% 
% %% Calculating the L2 norm of each image registration - registration
% 
% L2_norm_def = DEF_def;
% mean_Val = [];
% 
% 
% for i = 1:iLim
%     for j = 1:160
%         for k = 1:160
%             % Registration Deformation field minus Model Deformation field
%             L2_norm_def(i).img(j,k) = sqrt((regIm_def(:,i).img(j,k) - regIm_def(:,i).img(j,k))^2);
%             
%             if rDSMI_def(:,i).img(j,k) == 0
%                 L2_norm_def(i).img(j,k) = 0;
%             end      
%         end 
%     end
% end
% 
% %% Other statistical measures
% 
% % measure the average error value per image
% 
% mean_Val = [];
% std_Val = [];
% for i =1:iLim
%     
%     SimpleArray = L2_norm_def(:,i).img(:,:,1,1,2);
%     NewSimpleArray = SimpleArray(SimpleArray ~= 0);
%     mean_Val(i) = nanmean(NewSimpleArray,'all');
%     std_Val(i) = nanstd(NewSimpleArray);
%     
% end
% 
% %%
% 
% maskyy1 = rDSMI_def(:,1).img(:,:,1,1,1) - modDef(:,1).img(:,:,1,1,1);
% maskyy1(isnan(rDSMI_def(:,1).img(:,:,1,1,1) - modDef(:,1).img(:,:,1,1,1))) = 0;
% meany1 = mean(maskyy1,'all');
% 
% maskyy2 = rDSMI_def(:,1).img(:,:,1,1,2) - modDef(:,1).img(:,:,1,1,2);
% maskyy2(isnan(rDSMI_def(:,1).img(:,:,1,1,2) - modDef(:,1).img(:,:,1,1,2))) = 0;
% meany2 = mean(maskyy2,'all');
% 
% meanyLinear = mean(meany1+meany2);
% 
% 
% %% binary image
% 
% % mask2 = DEF_def2_mask(:,1).img(:,:,1,1,2);
% % % mask2(resamp_img2>0) = 1;
% % mask2(isnan(DEF_def2_mask(:,1).img(:,:,1,1,1))) = 0;
% 
% 
% 
% %%
% % meany = mean(mask2,'all');
% 
% 
% 
% %%
% stdVAl = std(mean_Val);
% 
% 
% % Measure the error in the Superior Inferior (SI) position:
% 
% 
% % Measure the error in the Anterior Posterior (AP) position
% %% Plot arbitrary figures to see if its working
% 
% figure(3);
% subplot(1,3,1)
% dispNiiSlice(L2_norm(:,12),"z",1,[-2 6])
% colorbar
% subplot(1,3,2)
% dispNiiSlice(L2_norm2(:,12),"z",1,[-2 6])
% colorbar
% subplot(1,3,3)
% dispNiiSlice(L2_norm3(:,12),"z",1,[-2 6])
% colorbar
% subplot(2,2,1)
% plot(mean_Val,x_20(101:115,1),'.');
% %% Plot Average value against corresponding signal
% t=1:iLim;
% figure(5);
% subplot(1,2,1)
% p = polyfit(t,std_Val,1);
% f = polyval(p,t);
% plot(t, std_Val,'.',t,f,'-');
% hold on 
% plot(t,f,'r-');
% % title('Mean pixel intensity value per image after registration against time')
% % xlabel('time (per image)');
% % ylabel('mean Pixel instensity value');
% % subplot(1,2,2)
% % plot(x_20(101:1500,1),mean_Val,'.');
% % xlabel('Surrogate signal');
% % ylabel('Pixel instensity value');
% 
% % 
% % subplot(1,3,3)
% plot(mean_Val,x_20(101:115,1),'.');