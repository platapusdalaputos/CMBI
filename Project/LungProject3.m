%% 6.5 Deformation Field Error (DFE)

% Pre-requisites

% Signed Distance Map Segmentation = SDMS
SDMS = load_untouch_nii('0007_sdt.nii.gz');

% Source mask
source_mask = load_untouch_nii('0007_mask.nii.gz');

% Source image
source_image = load_untouch_nii('0007.nii.gz');

%% Registration Images

% Deform the source image with the registration transformations
% Registation images = regIm
Z = 1;
for i = 101:1500
    
        [regIm(Z), regIm_def(Z), regIm_dis(Z)] = ...
            deformNiiWithCPGsSliding(cpg1(:,i), cpg2(:,i), SDMS, source_image, images(:,i));

        Z = Z+1;
end

%% Masking
% Deform the source mask image with the registration transformations
% The deformNiiWithCPGsSliding will also calculate the Deformation Field
% Signed Distance Map Segmentation = SDMS
% Source_image = SI
% Target_image = TI
% [The_Image , deformation_Field, distortion_Field] = ...
% ...deformNiiWithCPGsSliding(registration_Region1, registration_Region2, SDMS, SI, TI);

% Creating the registration Deformed Source Map Images (rDSMI)

Z = 1;
for i = 101:1500
    
        [rDSMI(Z), rDSMI_def(Z), rDSMI_dis(Z)] = ...
            deformNiiWithCPGsSliding(cpg1(:,i), cpg2(:,i), SDMS, source_mask, images(:,i));

        Z = Z+1;
end

% Creating the model Deformed Source Map Images (mDSMI)



% Multiplying the registration images by the registration...
% ...deformed mask source images creating: 
% Variable name given: New registration Masked images (NRMI)
%% Map the mask onto the registration images

% Create New registration Masked images (NRMI)

NRMI = regIm;

for i = 1:1400
    for j = 1:160
        for k = 1:160
            if rDSMI(:,i).img(j,k) == 0
                NRMI(:,i).img(j,k) = 0;
                rDSMI_def(:,i).img(j,k) = 0;                
            end      
        end 
    end
end


%% Plot arbitrary figures to see if its working


figure;
subplot(1,3,1)
dispNiiSlice(NRMI(:,108),"z",1)

subplot(1,3,2)
dispNiiSlice(regIm(:,108),"z",1)

subplot(1,3,3)
dispNiiSlice(rDSMI_def(:,108),"z",1)

%% Copy the header from all region 1 and region 2

Z = 1;
for i = 1:1400
    hdr1(Z) = cpg1(:,i).hdr;
    hdr2(Z) = cpg2(:,i).hdr;
    Z = Z+1;
end 

%% retrieve the signal for the 1400 images and input into linear model


surrogateSignals = [x_20(101:1500),ones(length(x_20(101:1500)),1)];
coefficients = C;
transformations = surrogateSignals*coefficients;



%% Split into four columns and reshape
cpg1_1 = reshape(transformations(1:1400,1:4489)', [67,67,1400,1,1]);
cpg1_2 = reshape(transformations(1:1400,4490:8978)', [67,67,1400,1,1]);
cpg2_1 = reshape(transformations(1:1400,8979:13467)', [67,67,1400,1,1]);
cpg2_2 = reshape(transformations(1:1400,13468:17956)', [67,67,1400,1,1]);

% cpg1_1 = transformations(1:1400,1:4489)';
% cpg1_2 = transformations(1:1400,4490:8978)';
% cpg2_1 = transformations(1:1400,8979:13467)';
% cpg2_2 = transformations(1:1400,13468:17956)';


%%
% cpg1_1a = cpg1_1;cpg1_2];
cpg1reshapefinal = zeros(67,67,1400,1,2);
cpg2reshapefinal = zeros(67,67,1400,1,2);

% 
for i = 1:1400
    cpg1reshapefinal(:,:,i,1,1)=cpg1_1(:,:,i);
    cpg1reshapefinal(:,:,i,1,2)=cpg1_2(:,:,i);
    cpg2reshapefinal(:,:,i,1,1)=cpg2_1(:,:,i);
    cpg2reshapefinal(:,:,i,1,2)=cpg2_2(:,:,i);
end

%% Create the transformations structs

Z = 1;
for i = 1:1400
    cpg1_new(Z).img = cpg1reshapefinal(:,:,i,:,:);
    cpg1_new(Z).hdr = hdr1(:,i);
    cpg2_new(Z).img = cpg2reshapefinal(:,:,i,:,:);
    cpg2_new(Z).hdr = hdr2(:,i);
    Z=Z+1;
end

%% Model Images

% Deform the source image with the Model transformations

Z = 1;
for i = 1:1400
    
        [modImg(Z), modDef(Z), modDis(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_new(:,i), cpg2_new(:,i), SDMS, source_image, images(:,i+100));

        Z = Z+1;
end
%% Masking
% Deform the source mask image with the model transformations
% Signed Distance Map Segmentation = SDMS
% Source_image = SI
% Target_image = TI
% [The_Image , deformation_Field, distortion_Field] = ...
% ...deformNiiWithCPGsSliding(registration_Region1, registration_Region2, SDMS, SI, TI);

% Creating the model Deformed Source Map Images (mDSMI)

Z = 1;
for i = i:1400
    
        [mDSMI(Z), mDSMI_def(Z), mDSMI_dis(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_new(:,i), cpg2_new(:,i), SDMS, source_mask, images(:,i+100));

        Z = Z+1;
end

%% Map the mask onto the modelled images

% Create New registration Masked images (NRMI)


for i = 1:1400
    for j = 1:160
        for k = 1:160
            if mDSMI(:,i).img(j,k) == 0
                modImg(:,i).img(j,k) = 0;
                modDef(:,i).img(j,k) = 0;                
            end      
        end 
    end
end

%% Plot arbitrary figures to see if its working


figure;
subplot(1,3,1)
dispNiiSlice(modImg(:,1200),"z",1)

% subplot(1,3,2)
% dispNiiSlice(mDSMI(:,1200),"z",1)

subplot(1,3,3)
dispNiiSlice(modDef(:,1200),"z",1)



    






