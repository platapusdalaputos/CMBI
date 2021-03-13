%% Functions for calculating the maximum tstat in each image pair

function MaxTstat = maxtstat(image1, image2, maskImage)

fid = fopen(image1, 'r', 'l'); % little-endian
data = fread(fid, 'float'); % 16-bit floating point
data = reshape(data, [40 40 40]); % dimension 40x40x40


fid1 = fopen(maskImage, 'r', 'l'); % little-endian
data1 = fread(fid1, 'float'); % 16-bit floating point
wm_mask = reshape(data1, [40 40 40]); % dimension 40x40x40

fid2 = fopen(image2, 'r', 'l'); % little-endian
data2 = fread(fid2, 'float'); % 16-bit floating point
data2 = reshape(data2, [40 40 40]); % dimension 40x40x40


CPA = data;
PPA = data2;

% Mapping the wm_mask onto the CPA images
% CPA mapping

for i=1:40
    for j=1:40
        for k=1:40 
             if wm_mask(k,j,i) == 0
                 CPA(k,j,i) = 0;
                                
            end
        end
    end    
end 

% Mapping the wm_mask onto the PPA images
% PPA mapping

for i=1:40
    for j=1:40
        for k=1:40 
             if wm_mask(k,j,i) == 0
                 PPA(k,j,i) = 0;
                                
            end
        end
    end    
end


% Computing the tstat on every voxel

tstat_value = [];
y1hodl = [];
y2hodl = [];
for i=1:40
    for j=1:40
        for k=1:40 
                 y1hodl(k) = CPA(k,j,i);
                 y2hodl(k) = PPA(k,j,i);

%              if k>39
                 transY1 = y1hodl';
                 transY2 = y2hodl';

                 tstat_value(k,j,i) = tstatt(transY1,transY2);
%                  if isnan(tstat_value(k,j,i))
%                      tstat_value(k,j,i)=0;
%                  end
%              end 
        end
    end    
end 


MakeSingleArray = tstat_value(tstat_value ~= 0);
MaxTstat = max(MakeSingleArray);




end

