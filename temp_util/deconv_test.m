clc;clear;
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]};
% neuron 4 134 136
addpath('E:\wupeixuan\CalImgProcess\CalciumPlot\utils');
load('E:\wupeixuan\auc_plot\data\Esr29m\Esr290201(no behave annnot)\Esr290201_8_combined.mat');
% trace0=trace(111,:);
% trace1=wdenoise(trace0,1); %���� trace0,�Ƿ����С�����룬foopsi��Ч������
trace1 = trace(109,:);
tic;
[c_cvx, s_cvx] = constrained_foopsi(trace1); 
toc;
tic
[c_cvx0, s_cvx0] = foopsi(trace1); 
figure('name', 'foopsi & constrained foopsi'); 
plot(trace1); hold on; plot(c_cvx); hold on; plot(c_cvx0);

%% AR1, foopsi
sprintf('AR1, foopsi:')
lambda = 1.5;
try
    tic;
    [dec_trace1, signal1, ~] = deconvolveCa(trace1, 'ar1', 'foopsi', 'lambda', lambda, 'optimize_pars');   % 0.130251
    toc;
catch
    fprintf('fail!\n');
end
try
    tic;
    [dec_trace1_0, signal1_0, ~] = deconvolveCa(trace1, 'ar1', 'foopsi', 'lambda', lambda); 
    toc;
    figure('name', 'AR1 foopsi');
    subplot(2,1,1);set(gca,'Position',[0.02,0.5,0.95,0.48]);
    plot(trace1); hold on; plot(dec_trace1);  % ��optimize����ϳ̶ȸ���
    title('AR1 foopsi with optimiz'); 
    subplot(2,1,2);set(gca,'Position',[0.02,0.01,0.95,0.48]);  % ��optimize�ĸ���exponential,�������˸���ϸ��
    plot(trace1); hold on; plot(dec_trace1_0);
    title('AR1 foopsi no optimiz');
catch
    fprintf('fail!\n');
end
%% AR2, foopsi
sprintf('AR2, foopsi:')
try
    tic;
    [dec_trace2, signal2, options] = deconvolveCa(trace1, 'ar2', 'foopsi', 'lambda', lambda, 'optimize_pars'); 
    toc;
    figure('name','AR2 foopsi');
    subplot(2,1,1);set(gca,'Position',[0.02,0.5,0.95,0.48]);
    plot(trace1); hold on; plot(dec_trace1_0);
    title('AR1 foopsi no optimiz'); 
    subplot(2,1,2);set(gca,'Position',[0.02,0.01,0.95,0.48]);  % ���߿���ȥû�����ϸ�ڱ����ԣ���Ҫ������
    plot(trace1); hold on; plot(dec_trace2);
    title('AR2 foopsi');
catch
    fprintf('fail!\n');
end
%% foopsi, conv kernal
sprintf('foopsi, conv kernal:')
g2 = [options.pars(1),  options.pars(2)];
taus = ar2exp(g2); 
w = 50;
ht = exp2kernel(taus, w); 
% % construct kernal with 2 exp % * �������
% sprintf('construct kernal with 2 exp:')
% lambda = 25; 
% try
%     tic;
%     [dec_trace3, signal3, options] = deconvolveCa(trace1, 'exp2', taus, 'foopsi', 'lambda', lambda, 'shift', w/2, 'window', w);
%     toc;
%     figure('name', 'FOOPSI, construct kernel');
%     plot(trace1); hold on; plot(dec_trace3);
% catch
%     fprintf('fail!\n');
% end
% use the kernel directly  % * trace4,5�����߶��뵽0��������
% sprintf('use the kernel directly:')
% try
%     tic;
%     [dec_trace4, signal4, ~] = deconvolveCa(trace1, 'kernel', ht, 'foopsi', 'lambda', lambda, 'shift', w/2, 'window', w);
%     toc;
% catch
%     fprintf('fail!\n');
% end
% % estimate the time constants
% sprintf('estimate the time constants:')
% lambda = 0; 
% try
%     tic;
%     [dec_trace5, signal5, ~] = deconvolveCa(trace1, 'exp2', 'foopsi', 'lambda', lambda, 'shift',  w/2, 'window', w, 'smin', 0.5);
%     toc;
%     figure('name', 'kernel & exp2');
%     subplot(2,1,1);set(gca,'Position',[0.02,0.5,0.95,0.48]);
%     plot(trace1); hold on; plot(dec_trace4);  % ��optimize����ϳ̶ȸ���
%     title('FOOPSI, kernel'); 
%     subplot(2,1,2);set(gca,'Position',[0.02,0.01,0.95,0.48]);  % ��optimize�ĸ���exponential
%     plot(trace1); hold on; plot(dec_trace5);
%     title('FOOPSI, exp2, estimate time const');
% catch
%     fprintf('fail!\n');
% end

%% constrained foopsi
sprintf('constrained foopsi:')
% estimate g with auto-correlation method
sprintf('estimate g with auto-correlation method:')
try
    tic;
    [dec_trace6, signal6, options] = deconvolveCa(trace1, 'ar1', 'constrained');   % ���Բ���ϸ��
    toc;
    g1 = options.pars;
    figure('name','Constrained_foopsi')
catch
    fprintf('fail!\n');
end
plot(trace1); hold on; plot(dec_trace6); hold on; plot(c_cvx);
% estimate g with auto-correlation method first and then update it to minimize the RSS
sprintf('estimate g with auto-correlation method first and then update it to minimize the RSS:')
try
    tic;
    [dec_trace7, signal7, ~] = deconvolveCa(trace1, 'ar1', 'constrained', 'optimize_b', 'optimize_pars');   % ! һ��Ҫ��'optimize_b'  ��dec_trace1������и��������\ϸ��, 0.577066s
    toc;
    figure('name','Constrained_foopsi, updating g')
    plot(trace1); hold on; plot(dec_trace7); hold on; plot(c_cvx);
catch
    fprintf('fail!\n');
end
%% AR1, hard thresh foopsi
sprintf('AR1, hard thresh foopsi:')
try
    tic;
    [dec_trace8, signal8, ~] = deconvolveCa(trace1, 'ar1', 'thresholded', 'optimize_smin', 'optimize_pars', 'thresh_factor', 0.99);  % �����˲���ϸ��
    toc;
    figure('name','AR1 thresh_foopsi with optimization')
    plot(trace1); hold on; plot(dec_trace8);
catch
    fprintf('fail!\n');
end
%% AR2, hard thresh foopsi
sprintf('AR2, hard thresh foopsi:')
try
    tic;
    [dec_trace9, signal9, ~] = deconvolveCa(trace1, 'ar2', 'thresholded','optimize_smin','optimize_pars', 'thresh_factor', 1);  % ���Բ���ϸ��
    toc;
    figure('name','AR2 thresh_foopsi with optimization');
    plot(trace1); hold on; plot(dec_trace9);
catch
    fprintf('fail!\n');
end
% %% AR2 kernel hard thresh foopsi  % ? �˴�����
% smin = 0.5; 
% sprintf('AR2 kernel hard thresh foopsi:')
% % exp
% sprintf('exp:')
% temp = roots([1, -g2(1), -g2(2)]);
% d = max(temp); 
% r = min(temp);
% pars = [d, r]; 
% try
%     tic;
%     [dec_trace10, signa10, ~] = deconvolveCa(trace1, 'exp2', pars, 'thresholded');
%     toc;
%     figure('name','AR2 exp kernel thresh_foopsi');
%     plot(trace1); hold on; plot(dec_trace10);
% catch
%     fprintf('fail!\n');
% end
% % kernel % * �ᱨ��
% sprintf('kernel:')
% ht = (exp(log(d)*(1:w)) - exp(log(r)*(1:w))) / (d-r); 
% try
%     tic;
%     [dec_trace11, signa11, ~] = deconvolveCa(trace1, 'kernel', ht, 'thresholded','smin', smin);
%     toc;
%     figure('name','AR2 kernel thresh_foopsi');
%     plot(trace1); hold on; plot(dec_trace11);
% catch
%     fprintf('fail!\n');
% end
% %% MCMC  % * û��
% sprintf('MCMC:')
% try
%     tic;
%     [dec_trace12, signa12, ~] = deconvolveCa(trace1, 'mcmc');
%     toc;
%     figure('name','AR1 MCMC');
%     plot(trace1); hold on; plot(dec_trace12);
% catch
%     fprintf('fail!\n');
% end
% %% AR2 MCMC % * û��
% sprintf('AR2 MCMC:')
% params.p=2;
% try
%     tic;
%     [dec_trace13, signa13, ~] = deconvolveCa(trace1, 'mcmc', params);
%     toc;
%     figure('name','AR2 MCMC');
%     plot(trace1); hold on; plot(dec_trace13);
% catch
%     fprintf('fail!\n');
% end

% various={'Time','original','c_cvx0','c_cvx','dec_trace6','dec_trace7','dec_trace1_0','dec_trace1','dec_trace2','dec_trace8','dec_trace9'};
% time0 = trace(1,:)';
% trace_table=table(time0,trace1',c_cvx0',c_cvx,dec_trace6,dec_trace7,dec_trace1_0,dec_trace1,dec_trace2,dec_trace8,dec_trace9,'VariableNames',various);
% writetable(trace_table, '\\liying.cibr.ac.cn\public\Code\auROC\wpx_code\traces_plot_with_events\Esr29\Esr290201(no behave annnot)\neuron136_deconv.csv');