function M = maxtstat(g1,g2,wm)

data = zeros(40,40,40,8);
data1 = zeros(40,40,40,8);

N = length(g1);

for i = 1:N
    
    fid = fopen(g1(i), 'r', 'l'); % little-endian
    H = fread(fid, 'float'); % 16-bit floating point
    data(:,:,:,i) = reshape(H, [40 40 40]); % dimension 40x40x40
    
    fid1 = fopen(g2(i), 'r', 'l'); % little-endian
    H1 = fread(fid1, 'float'); % 16-bit floating point
    data1(:,:,:,i) = reshape(H1, [40 40 40]); % dimension 40x40x40
    
end

fid2 = fopen(wm, 'r', 'l'); % little-endian
data2 = fread(fid2, 'float'); % 16-bit floating point
wm_mask = reshape(data2, [40 40 40]); % dimension 40x40x40




% CPA and PPA mapping

% Mapping the wm_mask onto the CPA and PPA images

CPA = data;
PPA = data1;

for l = 1:8
    for i=1:40
        for j=1:40
            for k=1:40 
                 if wm_mask(k,j,i) == 0

                     CPA(k,j,i,l) = 0;
                     PPA(k,j,i,l) = 0;
                     

                end
            end
        end    
    end 
end


% Computing the tstat on every voxel

tstat_value = [];
ttest_tstat = [];
y1hodl = [];
y2hodl = [];


for i=1:40
    for j=1:40
        for k=1:40 
            for l =1:8


                 y1hodl(l) = CPA(k,j,i,l);
                 y2hodl(l) = PPA(k,j,i,l);


                 transY1 = y1hodl';
                 transY2 = y2hodl';

                 tstat_value(k,j,i) = tstatt(transY1,transY2);

                
            end
        end
    end    
end 


M = max(tstat_value,[],'all');

end