% Displaying all images and make an animation of MR images for task 6.1

indxs=[1:10,400:410,800:810,1300:1310];


for k = 1:length(indxs)
  disp(k);
  
      ind=indxs(k);
    txt = strcat('test image : ',int2str(ind));

  % Runs only through first 100 images as its a little slow.
  % Problably because it has to handle many windows open and closes.
  
  % Keep the next line otherwise matlab cant save the figure.
  f = figure;
%   subplot(1,3,1)
  dispNiiSlice(modImg(ind),'z',1);
      t= text(-140,-190,txt);
    t.Color= [1,1,1];

  title('Linear Model Image')
  drawnow;
%   subplot(1,3,2)
%   dispNiiSlice(DEF_def2_mask2(k),'z',1,[-2 8]);
%       t= text(-140,-190,txt);
%     t.Color= [1,1,1];
% 
%   title('2nd Polynomial DFE')
%   subplot(1,3,3)
%   drawnow;
%   dispNiiSlice(DEF_def2_mask2(k),'z',1,[-2 8]);
%       t= text(-140,-190,txt);
%     t.Color= [1,1,1];
% 
%   title('3rd Polynomial DFE')
%   drawnow;
  % Stole this code from google.  
  frame = getframe(f); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 
  % Write to the GIF File 
  if k == 1 
      imwrite(imind,cm,'DFE.gif','gif', 'Loopcount',inf); 
  else 
      imwrite(imind,cm,'DFE.gif','gif','WriteMode','append'); 
  end 
  clf
  close % keep this 'close' otherwise it wont close any figure window and you
  % will end up with N number of overlapping windows.
end
