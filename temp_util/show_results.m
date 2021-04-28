load('E:\wupeixuan\auc_plot\data\Esr29m\Esr290218\Esr290218_deconv.mat');
load('E:\wupeixuan\auc_plot\data\Esr29m\Esr290218\Esr290218_7_prime_combined.mat');
y = trace(2,:)';
c_oasis = dec_data(1,:)';
s_oasis = sig_data(1,:)';
Fs = 10;
% clc;clear;
% % 4 134 136
% load('E:\wupeixuan\auc_plot\data\Esr29m\Esr290201(no behave annnot)\Esr290201_8_combined.mat');
% % trace0=trace(111,:);
% % trace1=wdenoise(trace0,1); %对于 trace0,是否加入小波降噪，foopsi的效果相似
% y = trace(109,:);
% tic;
% [c_cvx, s_cvx] = constrained_foopsi(y); 
% toc;
% tic
% [c_cvx0, s_cvx0] = foopsi(y); 
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]}; % colors
plot_cvx = false; 
figure( 'papersize', [15, 4]); 
set(gcf, 'color', 'w', ...
    'defaultAxesFontSize', 20, ...) 
    'defaultlinelinewidth',2, ...
    'position', [0, 0, 100*get(gcf, 'papersize')], ...
    'paperposition', [0, 0, get(gcf, 'papersize')])

% c
axes('position', [.05, .57, .95, .37]);
hold on;
plot(y, 'color', col{8}/255);
alpha(.7);
% plot(true_c, 'color', col{3}/255, 'linewidth', 1.5);
plot(c_oasis, '-.', 'color', col{5}/255);
if plot_cvx && exist('c_cvx', 'var')
    plot(c_cvx, '-.', 'color', col{7}/255);
end
axis tight;
last_num = ceil(length(dec_data)/Fs/100)*100;
xticks([0:75:last_num]);xlim([0 last_num]);
set(gca, 'xticklabel', []);
set(gca, 'ytick', 0:2);
ylabel('Fluor.');
box off;
if plot_cvx
    legend('Data', 'Truth', 'OASIS', 'CVX', 'location', 'northeast', 'orientation', 'horizontal');
else
    legend('Data', 'raw', 'OASIS', 'location', 'northeast', 'orientation', 'horizontal');
end
% s
axes('position', [.05, .18, .95, .37]);
hold on;
% plot(true_s, 'color', col{3}/255, 'linewidth', 1.5);
plot(s_oasis, '-.', 'color', col{5}/255);
if plot_cvx && exist('s_cvx', 'var')
    plot(s_cvx, '-.', 'color', col{7}/255);
end
axis tight;
last_num = ceil(length(dec_data)/Fs/100)*100;
xticks([0:75:last_num]);xlim([0 last_num]);
set(gca, 'xticklabel', get(gca, 'xtick')/30);
% set(gca, 'ytick', [0,1]);
xlabel('Time [s]');
ylabel('Activity.');