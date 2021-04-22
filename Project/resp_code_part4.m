%% Code for 3rd order Polynomial 

% Model Images poly 2

% Deform the source image with the Model transformations

Z = 1;
for i = 1:iLim
    
        [modImg3(Z), modDef3(Z), modDis3(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_p3(:,i), cpg2_p3(:,i), SDMS, source_image, images(:,i+iM));

        Z = Z+1;
end
%% Masking

% Deform the source mask image with the model transformations

% Variable name = Deformed Source Map Images (mDSMI)

Z = 1;
for i = 1:iLim
    
        [mDSMI3(Z), mDSMI_def3(Z), mDSMI_dis3(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_p3(:,i), cpg2_p3(:,i), SDMS, source_mask, images(:,i+iM));

        Z = Z+1;
end

%% Map the mask onto the deformed modelled images


for i = 1:iLim
    for j = 1:160
        for k = 1:160
            if mDSMI3(:,i).img(j,k) == 0
                modImg3(:,i).img(j,k) = 0;
                modDef3(:,i).img(j,k) = 0;                
            end      
        end 
    end
end

% %% Plot arbitrary figures to see if its working
% 
% 
% figure(2);
% subplot(1,3,1)
% dispNiiSlice(modImg3(:,1200),"z",1)
% title('Model Deformed Masked MRI Image')
% 
% subplot(1,3,2)
% dispNiiSlice(mDSMI3(:,1200),"z",1)
% title('Model Deformed Mask Image')
% 
% subplot(1,3,3)
% dispNiiSlice(modDef3(:,1200),"z",1)
% title('Model Deformed Mask Deformation Field Image')



% %% Deformation Field Error Again
% 
% % Create a new struct with the same values as rDSMI for the for loop
% DEF_def23 = rDSMI_def;
% 
% % Calculate the deformation field error for every voxel/pixel among the...
% % ...1400 target imgs.
% 
% for i = 1:1400
%     for j = 1:160
%         for k = 1:160
%             % Registration Deformation field minus Model Deformation field
%          DEF_def23(i).img(j,k) = regIm_def(:,i).img(j,k) - modDef3(:,i).img(j,k);  
%          
%         end 
%     end
% end

% %% Making the deformation with just the registration mask
% DEF_def2_mask3 = DEF_def23;
% 
% for i = 1:1400
%     for j = 1:160
%         for k = 1:160
%             if rDSMI_def(:,i).img(j,k) == 0
%                 DEF_def2_mask3(i).img(j,k) = 0;
%             end      
%         end 
%     end
% end
% % %% plotting deformatino field error
% figure;
% subplot(1,3,1)
% dispNiiSlice(DEF_def2_mask(:,1000),"z",1,[-2 3])
% colorbar;
% subplot(1,3,2)
% dispNiiSlice(DEF_def2_mask2(:,1000),"z",1,[-2 3])
% colorbar;
% subplot(1,3,3)
% dispNiiSlice(DEF_def2_mask3(:,1000),"z",1,[-2 3])
% colorbar;

% %% L2 norm
% L2_norm3 = DEF_def2_mask3;
% % mean_Val = [];
% 
% 
% for i = 1:1400
%     for j = 1:160
%         for k = 1:160
%             % Registration Deformation field minus Model Deformation field
%             L2_norm3(i).img(j,k) = sqrt((regIm_def(:,i).img(j,k) - modDef3(:,i).img(j,k))^2);
%             
%             if rDSMI_def(:,i).img(j,k) == 0
%                 L2_norm3(i).img(j,k) = 0;
%             end      
%         end 
%     end
% end

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
    Aa2(:,:,i) = modDef3(:,i).img(:,:,1,1,1);
    
end
%% Calculating the L2 norm of each image again BB



Bb1 = zeros(160,160,10);
Bb2 = zeros(160,160,10);
L2Bb = zeros(160,160,10);
for i = 1:iLim

    
    Bb1(:,:,i) = regIm_def(:,i).img(:,:,1,1,2);
    Bb2(:,:,i) = modDef3(:,i).img(:,:,1,1,2);
    
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


l2_mean_total_poly3 = mean(L2Aa+L2Bb,'all');


%%
figure(4);
subplot(1,3,1)
dispNiiSlice(L2_norm(:,9),"z",1,[-2 2])
colorbar;
subplot(1,3,2)
dispNiiSlice(L2_norm2(:,9),"z",1,[-2 2])
colorbar;
subplot(1,3,3)
dispNiiSlice(L2_norm3(:,9),"z",1,[-2 2])
colorbar;