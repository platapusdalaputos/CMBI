function [l2_mean_total_poly2] = secondPoly(iM, cpg1_mod, cpg2_mod, images, rDSMI, regIm_def)

% Signed Distance Map Segmentation = SDMS
SDMS = load_untouch_nii('0007_sdt.nii.gz');

% Source mask
source_mask = load_untouch_nii('0007_mask.nii.gz');

% Source image
source_image = load_untouch_nii('0007.nii.gz');

iTest = iM+1; % testing image start
iLim = 1500-iM;

% Model Images poly 2

% Deform the source image with the Model transformations

Z = 1;
for i = 1:iLim
    
        [modImg2(Z), modDef2(Z), modDis2(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_mod(:,i), cpg2_mod(:,i), SDMS, source_image, images(:,i+iM));

        Z = Z+1;
end



for i = 1:iLim
    for j = 1:160
        for k = 1:160
            if rDSMI(:,i).img(j,k) == 0
                modImg2(:,i).img(j,k) = 0;
                modDef2(:,i).img(j,k) = 0;                
            end      
        end 
    end
end


% Masking

% Deform the source mask image with the model transformations

% Variable name = Deformed Source Map Images (mDSMI)

% % % Z = 1;
% % % for i = 1:iLim
% % %     
% % %         [mDSMI2(Z), mDSMI_def2(Z), mDSMI_dis2(Z)] = ...
% % %             deformNiiWithCPGsSliding(cpg1(:,i), cpg2(:,i), SDMS, source_mask, images(:,i+iM));
% % % 
% % %         Z = Z+1;
% % % end
% % % 
% % % % Map the mask onto the deformed modelled images
% % % 
% % % 
% % % for i = 1:iLim
% % %     for j = 1:160
% % %         for k = 1:160
% % %             if mDSMI2(:,i).img(j,k) == 0
% % %                 modImg2(:,i).img(j,k) = 0;
% % %                 modDef2(:,i).img(j,k) = 0;                
% % %             end      
% % %         end 
% % %     end
% % % end


% Deformation Field Error Again

% Create a new struct with the same values as rDSMI for the for loop
% DEF_def22 = rDSMI_def;

% Calculate the deformation field error for every voxel/pixel among the...
% ...1400 target imgs.

% % % for i = 1:iLim
% % %     for j = 1:160
% % %         for k = 1:160
% % %             % Registration Deformation field minus Model Deformation field
% % %          DEF_def22(i).img(j,k) = regIm_def(:,i).img(j,k) - modDef2(:,i).img(j,k);  
% % %          
% % %         end 
% % %     end
% % % end

% Making the deformation with just the registration mask

% % % DEF_def2_mask2 = DEF_def22;
% % % 
% % % for i = 1:iLim
% % %     for j = 1:160
% % %         for k = 1:160
% % %             if rDSMI_def(:,i).img(j,k) == 0
% % %                 DEF_def2_mask2(i).img(j,k) = 0;
% % %             end      
% % %         end 
% % %     end
% % % end

% Calculate the l2 Norm of each image AA


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
% Calculate the L2 norm of each image BB



Bb1 = zeros(160,160,10);
Bb2 = zeros(160,160,10);
L2Bb = zeros(160,160,10);
for i = 1:iLim

    
    Bb1(:,:,i) = regIm_def(:,i).img(:,:,1,1,2);
    Bb2(:,:,i) = modDef2(:,i).img(:,:,1,1,2);
    
end

%

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

% Calculating the mean
% l2_mean_t = zeros(10,1);

for i = 1:iLim
    l2_mean_tAa(i,:) = mean(L2Aa(:,:,i),'all');
    l2_mean_tBb(i,:) = mean(L2Bb(:,:,i),'all');
end 


l2_mean_total_poly2 = mean(L2Aa+L2Bb,'all');

