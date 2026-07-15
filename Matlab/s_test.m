clc;clear;close all
path = 'D:\Code\Contraction_detection_manuscript_second_round\Contaction_detection\results\result\';
dir_list = dir(path);
dir_list = dir_list(3:end);
%%
f1_array = zeros(length(dir_list),12); %8 methods, 12 subjects
for i =1:length(dir_list)
    load([path dir_list(i).name]);
    f1_array(i,:) = f1_list;
end
show_stat(f1_array)
%%
group = {'G1','G2','G3','G4','G5','G6','G7','G8'};

figure(); % white background
boxplot(f1_array', group, ...           % black outlines
    'Symbol','o', ...               % outlier marker
    'Whisker',1.5, ...              % default whisker length
    'BoxStyle','outline', ...
    'Widths',0.6);

% Aesthetic tuning
set(gca,'FontSize',14,'LineWidth',1.2);
ylabel('F1 score','FontSize',16);
title('F1 score comparison','FontSize',16);

% Clean up look
box on;
set(gca,'GridLineStyle','--','GridAlpha',0.3);
%%
cp_score_array = zeros(length(dir_list),12); %8 methods, 12 subjects
for i =1:length(dir_list)
    load([path dir_list(i).name]);
    cp_score_array(i,:) = cp_score_list;
end
show_stat(cp_score_array)
%%
group = {'G1','G2','G3','G4','G5','G6','G7','G8'};

figure(); % white background
boxplot(cp_score_array', group, ...           % black outlines
    'Symbol','o', ...               % outlier marker
    'Whisker',1.5, ...              % default whisker length
    'BoxStyle','outline', ...
    'Widths',0.6);

% Aesthetic tuning
set(gca,'FontSize',14,'LineWidth',1.2);
ylabel('Scoring function error','FontSize',16);
title('Scoring function error comparison','FontSize',16);

% Clean up look
box on;
set(gca,'GridLineStyle','--','GridAlpha',0.3);
%%
cpn_array = zeros(length(dir_list),12); %8 methods, 12 subjects
for i =1:length(dir_list)
    load([path dir_list(i).name]);
    cpn_array(i,:) = cpn_list;
end
show_stat(cpn_array)
%%
group = {'G1','G2','G3','G4','G5','G6','G7','G8'};

figure(); % white background
boxplot(cpn_array', group, ...           % black outlines
    'Symbol','o', ...               % outlier marker
    'Whisker',1.5, ...              % default whisker length
    'BoxStyle','outline', ...
    'Widths',0.6);

% Aesthetic tuning
set(gca,'FontSize',14,'LineWidth',1.2);
ylabel('Change point number error','FontSize',16);
title('Change point number comparison','FontSize',16);

% Clean up look
box on;
set(gca,'GridLineStyle','--','GridAlpha',0.3);
%%
% x1 = f1_array(1,:);
% fprintf('median x1 = %.3f\n', median(x1));

%%
% for i = 1:length(dir_list)-1
%     x2 = f1_array(i+1,:);
%     fprintf('median x2 = %.3f\n', median(x2));
% %     [p, h, stats] = signrank(x1, x2,'alpha', 0.05, 'tail', 'right');
% %     
% %     % Display results
% %     fprintf('p-value = %.4f\n', p);
% %     fprintf('Test decision (h) = %d\n', h);
% %     disp(stats)
% % 
% % 
% %     disp('here')
%     alpha = 0.05;   % 95% CI
%     nboot = 10000;
%     [p,h,stats] = signrank(x1,x2);
%     
%     diffs = x1 - x2;
% %     HL = median(diffs);
% %     bootstat = bootstrp(nboot, @(d) median(d), diffs);
% %     ci = prctile(bootstat, [100*alpha/2, 100*(1-alpha/2)]);
%     diffs(diffs == 0) = [];
%     N = numel(diffs);
%     
%     % Compute Z statistic
%     W = stats.signedrank;
%     Z = (W - (N*(N+1))/4) / sqrt(N*(N+1)*(2*N+1)/24);
%     
%     % Compute effect size
%     r = Z / sqrt(N);
%     
%     fprintf('p = %.4f, Z = %.4f, r = %.4f\n', p, Z, abs(r));
% end
% 
% 
% 
% 
% 
