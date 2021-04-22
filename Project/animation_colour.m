% Displaying all images and make an animation of MR images for task 6.1

% load the source image and signe distance map
% source_nii = load_untouch_nii('0007.nii.gz');
% dist_nii = load_untouch_nii('0007_sdt.nii.gz');

% indices of images to load
indxs=(1:30);

for k = 1:length(indxs)
    disp(k);

    % Runs only through first 100 images as its a little slow.
    % Problably because it has to handle many windows open and closes.

    ind=indxs(k);

%     % deform an image with a B-spline transformation 
%     % Linear
%     [def_vol_nii_lin, def_field_nii_lin, dis_field_nii_lin] = ...
%     deformNiiWithCPGsSliding(cpg1_lin(ind), cpg2_lin(ind), dist_nii, source_nii, images(100+ind));
%     % 2nd Polynomial
%     [def_vol_nii_p2, def_field_nii_p2, dis_field_nii_p2] = ...
%     deformNiiWithCPGsSliding(cpg1_p2(ind), cpg2_p2(ind), dist_nii, source_nii, images(100+ind));
%     % 3rd Polynomial
%     [def_vol_nii_p3, def_field_nii_p3, dis_field_nii_p3] = ...
%     deformNiiWithCPGsSliding(cpg1_p3(ind), cpg2_p3(ind), dist_nii, source_nii, images(100+ind));

    % prepare text do display on figures
    txt = strcat('test image : ',int2str(ind));
        
    % Keep the next line otherwise matlab cant save the figure.
    f = figure;
    
    % plot linear
    
    subplot(2,2,1);
    dispNiiSliceColourOverlay(L2_norm(ind),L2_norm_def(ind),"z",1, [-2 8]); 
    t= text(-140,-190,txt);
    t.Color= [1,1,1];
    title("Lin mod vs Orgin");
    
    % plot p2
    subplot(2,2,2);
    dispNiiSliceColourOverlay(L2_norm2(ind),L2_norm_def(ind),"z",1, [-2 8]); 
    t= text(-140,-190,txt);
    t.Color= [1,1,1];
    title("2nd Poly. vs Orgin");
    
    % plot p3
    subplot(2,2,3);
    dispNiiSliceColourOverlay(L2_norm3(ind),L2_norm_def(ind),"z",1, [-2 8]);
    t= text(-140,-190,txt);
    t.Color= [1,1,1];
    title("3rd Poly. vs Orgin");
    
    % plot p4
    subplot(2,2,4);
    dispNiiSlice(modImg(ind),'z',1);
    t= text(-140,-190,txt);
    t.Color= [1,1,1];
    title('Linear Model Image');
  
    
    drawnow;

    % Stole this code from google.  
    frame = getframe(f); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if k == 1 
        imwrite(imind,cm,'modvsorig.gif','gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,'modvsorig.gif','gif','WriteMode','append'); 
    end 
    clf
    close % keep this 'close' otherwise it wont close any figure window and you
    % will end up with N number of overlapping windows.
end