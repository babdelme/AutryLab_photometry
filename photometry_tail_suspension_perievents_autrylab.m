%% Script for visualization and analysis of photometry traces during tail suspension 
%% Author: Ilaria Carta
% Photometry traces were acquired with Bonsai 2.4. The Bonsai workflow
% splits the original data (with interleaved GCamp and isosbestic datapoints) and outputs 2 traces.
% Behavior data was analyzed with Ethovision and the resultscexported as binary files.
% This code plots the 2 traces, corrects the data for bleaching, uses the isosbestic trace to normalize the data, computes deltaF/F, aligns corrected photomtrey data with behavior
% Finally, it plots percent change in deltaF/F and z-score.

clear all
close all


%% open 470 trace
filename = 'C:\Users\username\Documents\filename_470';              %extension not specified because Bonsai output files do not have extension
delimiter = ' ';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
data470 = dataArray{:, 1};
fclose(fileID);
clearvars filename delimiter formatSpec fileID dataArray ans;

%% open 410 trace

filename = 'C:\Users\username\Documents\filename_410';              %extension not specified because Bonsai output files do not have extension
delimiter = ' ';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
data410 = dataArray{:, 1};
fclose(fileID);
clearvars filename delimiter formatSpec fileID dataArray ans;
photometry_framerate = 20;
recording_duration_s = floor(length(data470)/photometry_framerate);


%% open behavior analysis file and specify parameters

behav_scoring = load('mouseID_behavior_filename.txt');
behavior_frame_rate=30;
video_duration_s = floor(length(behav_scoring)/behavior_frame_rate);
data_offset = video_duration_s - recording_duration_s;
behav_cut=behav_scoring(data_offset*behavior_frame_rate:end,:);

tail_suspension=behav_cut(:,4);

tailsusp_events = tail_suspension';



%% Define time and time bins

d=recording_duration_s;
time=1:length(data470);
start_point = 1;               % Enter the desired starting point in seconds, if you wish to crop the first part of the trace. 
first_bin=time(start_point*photometry_framerate:photometry_framerate*d);   %in this case it's the first and only time bin, but it is possible to define more time bins
data1=data470(start_point*photometry_framerate:photometry_framerate*d,1);
data2=data410(start_point*photometry_framerate:photometry_framerate*d,1);


%% plot photometry data only (1st figure)
figure
subplot (3,1,1);
plot(first_bin,data1, 'r');
xlabel('time(s)')
ylabel('avg pixel int')
xticks(1*photometry_framerate:startpoint*photometry_framerate:photometry_framerate*d)
xticklabels(0:30:d)
axis tight
ylim([78 81.5]);

subplot (3,1,2);
plot(first_bin, data2, 'Color', [0.5 0 0.5]);
xlabel('time(s)')
ylabel('avg pixel int')
xticks(1*photometry_framerate:startpoint*photometry_framerate:photometry_framerate*d)
xticklabels(0:30:d)
axis tight
ylim([6.7 7.7]);


%% Correct for bleaching, scale, and normalize the trace (this part was kindly provided by Sage Aronson, Neurophotometrics, and minimally adapted)
figure
subplot(2,1,1)
time_F_matrix=[first_bin' data1 data2];
temp_fit = fit(time_F_matrix(:,1),time_F_matrix(:,2),'exp2');
plot(first_bin,data1-temp_fit(first_bin));
title('Linearize by Fitting with Biexponential')
xlabel('Time of day in total ms')
ylabel('F')
clear temp_fit

subplot (2,1,2)
%fit isosbestic
temp_fit = fit(time_F_matrix(:,1),time_F_matrix(:,3),'exp2'); %note, in this case, I am ignoring the first 2000 points where there is this weird fast decay to get a better fit. experimentally, i normally set things up so this isn't an important time in the recording / animal is habituating.
%scale fit to calcium dependent data
fit2 = fitlm(temp_fit(time_F_matrix(:,1)),time_F_matrix(:,2));
%calculate a crude dF/F by subtracting and then dividing by the fit
plot(time_F_matrix(:,1),100*(time_F_matrix(:,2)-(fit2.Fitted))./(fit2.Fitted))
xlabel('Time(s)')
ylabel('crude dF/F (%)')
title('Linearizing + Normalizing Using Isosbestic')

%data correction
FP.fakebackground = 0.017;                                                  %fakebackround because it's only an estimate
FP.corrected_data(:,1) = time_F_matrix(:,1);                                %first column is the time vector

figure
for i = 2:3
       temp_fit = fit(time_F_matrix(:,1),time_F_matrix(:,3),'exp2'); 
%scale fit to calcium dependent data
    fit2 = fitlm(temp_fit(time_F_matrix(:,1)),time_F_matrix(:,i));
%calculate a crude dF/F by subtracting and then dividing by the fit
subplot(3,1,i);

%calculate dF/F by subtracting the background
    FP.corrected_data(:,i) = (time_F_matrix(:,i)-((fit2.Fitted)-FP.fakebackground))./(fit2.Fitted-FP.fakebackground); 
clear temp_fit fit2
end

f5 = figure;
for i = 1:2
    subplot(2,1,i)
   
    plot(FP.corrected_data(:,1),FP.corrected_data(:,i+1))
     title('Corrected Data')
     xlabel('time of day in total ms')
     ylabel('%dF/F')
end

GCamP_trace_corrected_smooth = smooth(FP.corrected_data(:,2),0.001);
GCamP_trace_abs = GCamP_trace_corrected_smooth - min(GCamP_trace_corrected_smooth);

%% plot data from photometry and behavior together (2nd figure)

figure
ax1=axes('units','inches','position',[2 2.6 12 4]);                        %in order to overlay two plots with different number of data points such as photometry and behavioral data, we need to specify the location of the 1st and second plot, this line specifies the first.
plot(first_bin,GCamP_trace_abs,'k','LineWidth',1.5);                                         
xticks(1*photometry_framerate:60*photometry_framerate:photometry_framerate*d);
xticklabels(0:60:d);
axis tight;
box off
ylabel('\DeltaF/F');
xlabel('time(s)');
ylim([-0.007 0.02]);
% xlim([100*photometry_framerate 500*photometry_framerate]);


%% plot behaviors underneath trace

ax2=axes('units','inches','position',[2 2.6 12 .2]);                         %this specifies the position of the second plot (behavior data exported from ethovision) 
imagesc(tailsusp_events(start_point*behavior_frame_rate:d*behavior_frame_rate),'AlphaData', 0.6);
colormap(ax2,[0.8 0.8 0.8; 0 1 0; 0 0 1]);  
set(gca,'YTick',[]);
set(gca,'XTick',[]);
box off
hold on


%% create ROIs with user input
% The following loop allows users to manually draw ROIs and saves variables within the chosen ROI

roi_s=[];                                                                  % ROI start
roi_e=[];                                                                  % ROI end
roi_d=[];                                                                  % ROI duration

a = input('Would you like to draw an ROI?\n(1/0)\n ');
fps = photometry_framerate;
j=1;
    while a==1
    fprintf ('please, draw a rectangle\n');

roi = drawrectangle(ax1);
pause
r=roi.Position;

roi_s(1,j) = r(1,1);
roi_e(1,j) = r(1,1)+r(1,3);
roi_d(1,j) = r(1,3);

px = [0 1 1 0]*roi_d(1,j)+roi_s(1,j);
py = [0 0 1 1]*1.5+(-0.15);
if j<2
patch(px, py, 'Red', 'FaceColor', [0 0 0.8], 'LineWidth', 0.1, 'LineStyle', '-', 'FaceAlpha', 0.2,'edgealpha',0.8);
else 
    patch(px, py, 'Red', 'FaceColor', [0.8 0 0.1], 'LineWidth', 0.1, 'LineStyle', '-', 'FaceAlpha', 0.2,'edgealpha',0.8);
end
qx = [0 1 1 0]*20+roi_e(1,j);
qy = [0 0 1 1]*1.5+(-0.15);
patch(qx, qy, 'Red', 'FaceColor', [0.8 0 0.1], 'LineWidth', 0.1, 'LineStyle', '-', 'FaceAlpha', 0.2,'edgealpha',0.8);

j = j+1;
    
fprintf ('ROI borders saved!')

a = input('Would you like to draw another ROI?\n(1/0)\n \n');
    end
if a==0
    fprintf ('All ROIs defined. Proceeding with analysis...');
else
    return
end

%ROI CALCULATIONS
%% this part performs some basic calculations (%change in deltaF/F and z-score) on the ROI data saved previously

baseline = GCamP_trace_abs(round(roi_s(1,1)):round(roi_s(1,1))+15*fps);    %set a 15s baseline
baseline_events = cell(length(roi_s),1);                                   %preallocate cell for next loop
events = cell(length(roi_s),1);
peri_events = cell(length(roi_s),1);

for i=1:length(roi_s)    
baseline_events{i,1} = [baseline' GCamP_trace_abs(round(roi_s(1,i)):round(roi_e(1,i)))'];  % segments containing baseline and event for each iteration. roi_1(1,1) is a different baseline, will be removed after
events{i,1} = GCamP_trace_abs(round(roi_s(1,i)):round(roi_e(1,i)));                        % includes first event(baseline), will be removed after
peri_events{i,1} = GCamP_trace_abs(round(roi_s(1,i)-10*fps):round(roi_e(1,i)))';           % considers the event plus the prior 10s as a whole segment, includes first event(baseline), will be removed after
end

all_events = peri_events{2:length(peri_events)};

% this part converts the event matrix into a cell array, replaces empty values with NaNs, then converts back into a matrix (adapted from Mathworks solutions)

maxSize = max(cellfun(@numel, peri_events));                               % Get the maximum vector size
fcn = @(x) [x nan(1, maxSize-numel(x))];                                   % Create a function to fill empty values (NaN)
rmat = cellfun(fcn, peri_events, 'UniformOutput', false);                  % Pad each cell with NaNs
rmat = vertcat(rmat{:});                                                   % arrange data for plotting the heatmap later 
events_stacked_matrix = rmat;
events_only = events_stacked_matrix(1:length(events), length(baseline):end)';
events_means_mouseID = nanmean(events_only);
events_means_mouseID_except_baseline = events_means_mouseID(2:end);        % the manually drawn baseline is not used if we are calculating perievents, but may be useful for long lasting changes (predator odor, forced swim test, etc.)
avg_events_except_baseline = mean(events_means_mouseID_except_baseline);

means_to_compare_mouseID_events = [events_means_mouseID(1,1), avg_events_except_baseline];    %will be used later for plotting

% same for the perievent matrix
all_peri_events = peri_events{2:length(peri_events)};

maxSize = max(cellfun(@numel, peri_events));                               % Get the maximum vector size
fcn2 = @(x) [x nan(1, maxSize-numel(x))];                                  % Create a function to fill empty values (NaN)
rmat2 = cellfun(fcn, peri_events, 'UniformOutput', false);                 % Pad each cell with NaNs
rmat2 = vertcat(rmat2{:});                                                 % arrange data for plotting the heatmap later 
peri_events_stacked_matrix = rmat2;
peri_events_only = events_stacked_matrix(2:length(events),:)';
pre_events_only = peri_events_only(1:10*fps,1:length(events)-1);
pre_events_means_mouseID = nanmean(pre_events_only);
post_events_only = peri_events_only(10*fps+1:end,1:length(events)-1);
post_events_means_mouseID = nanmean(post_events_only);
avg_pre_events = mean(pre_events_means_mouseID);
avg_post_events = mean(post_events_means_mouseID);
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));  %function created to calculate zscores of a matrix ignoring empty values (NaN)

z_mouseID_events = zeros(length(roi_s),length(events_stacked_matrix));     %will be used later for plotting
z_mouseID_peri_events = zeros(length(roi_s)-1,length(peri_events_stacked_matrix));  %will be used later for plotting

means_to_compare_mouseID_peri_events = [avg_pre_events, avg_post_events]; %will be used later for plotting

for i = 1:length(roi_s)
    z_mouseID_events(i,:) = zscor_xnan(events_stacked_matrix(i,:));        %calculates z-scores of the entire matrix  
end

for i = 1:length(roi_s)
    z_mouseID_peri_events(i,:) = zscor_xnan(peri_events_stacked_matrix(i,:));        %calculates z-scores of the entire matrix  
end
%% the following figure plots percent change in deltaF/F and z-score

figure
subplot (3,2,1)
bar((means_to_compare_mouseID_events/events_means_mouseID(1,1))*100);      %mean of all events against mean of the manually drawn baseline. here we take the manually drawn baseline, but can be adapted to use the pre-event baseline
ylabel('% \DeltaF/F');
ylim([0 200]);
hold on
plot(xlim,[100 100],'LineStyle','--','Color','r','LineWidth',1.5);
set(gca,'XTick',[]);
box off

subplot (3,2,2)
bar((means_to_compare_mouseID_peri_events/avg_pre_events)*100);            %mean of all post events against mean of all pre-events
ylim([0 200]);
ylabel('% \DeltaF/F');
hold on
plot(xlim,[100 100],'LineStyle','--','Color','r','LineWidth',1.5);
set(gca,'XTick',[]);
box off

subplot (3,2,3)
bar((events_means_mouseID/events_means_mouseID(1,1))*100);                 %mean of each individual event against mean of the manually drawn baseline. 
ylim([0 700]);
ylabel('% \DeltaF/F');
hold on
plot(xlim,[100 100],'LineStyle','--','Color','r','LineWidth',1.5);
set(gca,'XTick',[]);
box off

subplot (3,2,4)                                                            %mean of each individual post events against mean of each individual pre-events
pre_post_means_to_compare = [avg_pre_events post_events_means_mouseID];
bar((pre_post_means_to_compare/avg_pre_events)*100);
ylim([0 200]);
ylabel('% \DeltaF/F');
hold on
plot(xlim,[100 100],'LineStyle','--','Color','r','LineWidth',1.5);
set(gca,'XTick',[]);
box off 

subplot (3,2,5)                                                            %displays heatmap of z-scores where each row is an event, each column the zscore of a datapoint in time.
imagesc(z_mouseID_events(2:end,:))
cb=colorbar;
xticks(1:5*fps:20*fps);
xticklabels(-10:5:10);
xlim([0 400]);
ylabel('#mouseID');
xlabel('time from event onset(s)');
set(get(cb,'Title'),'String','z-score');   
box off 

subplot (3,2,6)
imagesc(z_mouseID_peri_events(2:end,:))                                   %displays heatmap of z-scores where each row is a perievent, each column the zscore of a datapoint in time.
cb=colorbar;
xticks(1:10*fps:100*fps);
xticklabels(-10:10:100);
ylabel('#mouseID');
xlabel('time from event onset(s)');
set(get(cb,'Title'),'String','z-score');   
box off
