%% Script for averaging data generated using "photometry_tail_suspension_perievents_autrylab"
%% Author: Ilaria Carta
% This code works with previously saved variables. The variables were generated with the script
% "photometry_tail_suspension_perievents_autrylab" and saved in a .mat file.
% Given those variables, this code computes and plots average zscores, average percent change deltaF/F and area under the curve.

clear all
close all
load('GalCre_mice_jan22_female_exp.mat');

%% calculate zscores using baseline zscore
peri_event_variables = who('-regexp','peri_events_stacked_matrix_*');      %finds files containing the specified characters in quotes
N=numel(peri_event_variables);                                             %counts the number of files containing the specified characters, for averaging purposes
length_peri_event = 600;                                                   %this and the following two lines allow to specify the desired perievent duration in s
length_baseline = 200;
length_event = 400;
segments=nan(N,600);                                                       %preallocate variables before loop
baseline_10s = nan(N,length_baseline);                                     %             ||
zscore_container=nan(N,length_peri_event);                                 %             ||

%% for loop to crop events to the same dimension and compute z-score
for i=1:N
    e=eval(peri_event_variables{i,1});
    segments(i,:)=e(2,1:length_peri_event);
    baseline_10s(i,:) = segments(i,1:length_baseline);
    zscore_container(i,:)=(segments(i,:)-mean(baseline_10s(i,:)))./std(baseline_10s(i,:));    
end

%% plot data
fps=20;                                                                     %specify frame rate

figure

subplot (3,2,1)
im_data = imagesc(zscore_container);                                        %plots zscore in the form of a heatmap preserving all individual mouse data (rows)
h = im_data;
h.AlphaData = ones(size(h.CData));
h.AlphaData(isnan(h.CData)) = 0;                                            %NaNs will show white in the heatmap (optional)
cb=colorbar;
caxis([0 20]);                                                              %caxis decides the lower and upper limits of the colorbar. 
xticks(1:5*fps:50*fps);
xticklabels(-10:5:20);
row_titles = [  'mouseID';  'mouseID';  'mouseID';  'mouseID'];             %keep the same number of characters for each animal
yticks(1:N);                                                               
yticklabels({row_titles});
xlabel('time from approach(s)');
set(get(cb,'Title'),'String','z-score'); 
xlim([1 600])

subplot (3,2,5)

means_stacked = [means_to_compare_mouseID_events;means_to_compare_mouseID_events;means_to_compare_mouseID_t1_events;means_to_compare_mouseID_t1_events];  %specify mouseID
normalized_avg=[(means_stacked(:,1)./means_stacked(:,1))*100, (means_stacked(:,2)./means_stacked(:,1))*100];
means_all = mean(normalized_avg);

bar(means_all);                                                            %plots average change in deltaF\F          
ylim([0 300]);
ylabel('% \DeltaF/F');
hold on
plot(xlim,[100 100],'LineStyle','--','Color','r','LineWidth',1.5);
hold on
plot(1:2,normalized_avg,'o');
xlab = {'pre';'post'};
set(gca,'xtick',1:2,'xticklabel',xlab, 'Fontsize' , 8);
box off

subplot (3,2,3)

alpha = 0.15;
stdshade(zscore_container, alpha, 'b');                                    %plots average zscore with standard deviation (requires the stdshade function, available on Matworks)
axis tight
xticks(1:5*fps:50*fps);
xticklabels(-10:5:20); 
xlabel('time from approach(s)');
ylabel('z-score');
title('avg z score') 
xlim([1 600])
ylim([-10 50])
%ylim([-2 2.5]);

   
%% compute area under the curve with the built-in function trapz

rawdata =[];
rawdata(1,:)=peri_events_stacked_matrix_mouseID_t1(2,1:600);               %specify mouseID
rawdata(2,:)=peri_events_stacked_matrix_mouseID_t1(2,1:600);               %      ||
rawdata(3,:)=peri_events_stacked_matrix_mouseID_t1(2,1:600);               %      ||
rawdata(4,:)=peri_events_stacked_matrix_mouseID_t1(2,1:600);               %      ||
c_1st_half_1 = (rawdata(1, 1:200))+abs(min(rawdata(1, 1:200)));t1_1=mean(c_1st_half_1);
c_2nd_half_1 = (rawdata(1, 201:400))+abs(min(rawdata(1, 1:200)));t1_2=mean(c_2nd_half_1);
c_1st_half_2 = (rawdata(2, 1:200))+abs(min(rawdata(2, 1:200)));t2_1=trapz(c_1st_half_2);
c_2nd_half_2 = (rawdata(2, 201:400))+abs(min(rawdata(2, 1:200)));t2_2=trapz(c_2nd_half_2);
c_1st_half_3 = (rawdata(3, 1:200))+abs(min(rawdata(3, 1:200)));t3_1=trapz(c_1st_half_3);
c_2nd_half_3 = (rawdata(3, 201:400))+abs(min(rawdata(3, 1:200)));t3_2=trapz(c_2nd_half_3);
c_1st_half_4 = (rawdata(4, 1:200))+abs(min(rawdata(4, 1:200)));t4_1=trapz(c_1st_half_4);
c_2nd_half_4 = (rawdata(4, 201:400))+abs(min(rawdata(4, 1:200)));t4_2=trapz(c_2nd_half_4);


T_pre = [t1_1 t2_1 t3_1 t4_1];
T_post = [t1_2 t2_2 t3_2 t4_2];

subplot (3,2,6)                                                            % plots average AUC
auc_rawdata = [T_pre' T_post'];
bar(1:2,mean(auc_rawdata));
hold on
plot(1:2,auc_rawdata,'o');
xlab = {'pre';'post'};
set(gca,'xtick',1:2,'xticklabel',xlab, 'Fontsize' , 8);
ylabel('AUC');
box off

