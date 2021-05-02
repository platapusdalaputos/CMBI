function [inhalation, exhalation] = inh_exh(n)

% Splitting the surrogate signals

inhalation = [];
exhalation = [];

for i = 2:n
    if x_20(i) > x_20(i-1)
        inhalation(i) = x_value(i);
%         inhalation(i+1) = x_20(i);
    else
        exhalation(i) = x_value(i);
%         exhalation(i+1) = x_20(i);
    end
end
   

inhale = nonzeros(inhalation)';
inhalation(inhalation==0) = NaN;
exhale = nonzeros(exhalation)';
exhalation(exhalation==0) = NaN;


% Plot the figure for inhale_exhale 

figure;

plot(inhalation(1:75),'-o','MarkerFaceColor', 'b','LineWidth',2)
hold on
plot(exhalation(1:75),'-o','MarkerFaceColor', 'r','LineWidth',2)
ax = gca;
xlabel('Frame','FontSize', 16);
ylabel('Pixel','FontSize', 16);
legend('Inhalation', 'Exhalation','Location','northwest')
        
exportgraphics(ax,fullfile('./figures', sprintf('%s_inhale_exhale.png',x_value)))


