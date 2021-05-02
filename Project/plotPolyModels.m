function [] = plotPolyModels(iM, x_value1, coeff, coeff2, coeff3, Sig2, Sig3, SI_deform_CP_44_38, S_lin, title)

% retrieve the coefficient of CP Region 1 SI and reshape them in a 67x67 grid
C_reg1_SI = reshape(coeff(:,(4489+1):8978),[],67,67);
C_p2_reg1_SI = reshape(coeff2(:,(4489+1):8978),[],67,67);
C_p3_reg1_SI = reshape(coeff3(:,(4489+1):8978),[],67,67);

% get the corresponding coefficients
C_44_38 = C_reg1_SI(:,44,38);
C_p2_44_38 = C_p2_reg1_SI(:,44,38);
C_p3_44_38 = C_p3_reg1_SI(:,44,38);


% plot the estimated models at CP 44-38

figure;
plot(x_value1(1:iM),SI_deform_CP_44_38(1:iM),':')
hold on;
plot(x_value1(1:iM),SI_deform_CP_44_38(1:iM),'bx');
hold on;
plot(x_value1(1:iM),S_lin*C_44_38,'r-');
hold on
[sorted_p2,ind] = sort(Sig2*C_p2_44_38,'descend');
plot(x_value1(ind), sorted_p2,'g-');
hold on
[sorted_p3,ind] = sort(Sig3*C_p3_44_38,'descend');
plot(x_value1(ind), sorted_p3,'m-');
ax = gca; 
ylim([30,58])        
xlabel("surrogate value",'FontSize', 16);
ylabel("control-point value",'FontSize', 16)
legend(' ',' ','linear', 'Poly^{2nd}', 'Poly^{3rd}');
% legend(' ',' ',

exportgraphics(ax,fullfile('./figures', sprintf('%s_polymodel_%d_trainingImages.png',title,iM)))


