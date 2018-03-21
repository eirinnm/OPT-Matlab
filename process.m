videoname = 'Mut-2017-09-08T19_15_12.avi';
video = read(VideoReader(videoname));
background = imread('background.bmp');
video = (background - video);
disp('Video loaded');
%%
% grab a frame
ff = video(:,:,1,20);
top = sgolayfilt(double(mean(ff(1:10,:))),5,55);
%bottom = -single(mean(ff(end-10:end,:)));
findpeaks(top, 'MinPeakProminence', 10, 'MinPeakDistance',300);
%[pks, loc] = findpeaks(bottom, 'MinPeakProminence', 40, 'SortStr','descend');
%pks=pks(1:2);
%loc=loc(1:2);
%% Get the sinogram of the first row
%topsino = -single(squeeze(video(1,:,1,:)));
topsino = double(squeeze(mean(video(1:10,:,1,:),1)));
bottomsino = double(squeeze(mean(video(end-10:end,:,1,:),1)));
topmid = [];
toppoints = [];
bottompoints = [];
bottommid = [];
for t = 1:size(video,4)
    top_smoothed = sgolayfilt(topsino(:,t),5,55);
    [pks,loc] = findpeaks(top_smoothed, 'MinPeakProminence', 10,'MinPeakDistance',300, 'SortStr','descend');
    % get the first two peaks by height and average their position
    if size(loc,1)>=2
        midpoint = mean(loc(1:2));
        toppoints = [toppoints; loc(1:2)'];
    else
        midpoint = nan;
        toppoints = [toppoints; [nan, nan]];
    end
    topmid = [topmid; midpoint];
    %now do the same thing for the bottom sinogram
    bottom_smoothed = sgolayfilt(double(bottomsino(:,t)),5,55);
    [pks,loc] = findpeaks(bottom_smoothed, 'MinPeakProminence', 10,'MinPeakDistance',300, 'SortStr','descend');
    if size(loc,1)>=2
        midpoint = mean(loc(1:2));
        bottompoints = [bottompoints; loc(1:2)'];
    else
        midpoint = nan;
        bottompoints = [bottompoints; [nan, nan]];
    end
    bottommid = [bottommid; midpoint];
end
imagesc(bottomsino);
set(gca,'YDir','normal');
hold on;
plot(bottompoints);
hold off;

%%
%replace error points with nans
%figure;
topmid(abs(diff(topmid))>20) = nan;
bottommid(abs(diff(bottommid))>20) = nan;
plot(topmid);
hold on;
plot(bottommid);
hold off;
title('Midpoint of top and bottom capillary walls');
xlabel('Frame number');
ylabel('Pixel');
%% Get the COR from these points
midpoints = [topmid, bottommid];
COR = nanmean(midpoints,1);
%COR = [nanmean(topmid), nanmean(bottommid)];
figure;
ff = video(:,:,1,125);
imagesc(ff);
hold on;
plot(COR,[1, size(video,1)]);
hold off;
%% How many frames for a full rotation?
startframe_return = [];
for startframe = 1:20
    d = abs(midpoints(startframe,:) - midpoints(300:end,:));
    [m, endframe] = min(sum(d,2));
    startframe_return = [startframe_return; [m, endframe]];
end
[m,best_startframe]=min(startframe_return(:,1));
best_endframe = 300 + startframe_return(best_startframe,2);
% generate sequence of angles - assuming it's uniform
angles = linspace(0,2*pi, best_endframe-best_startframe+1);

%% Clip the video and centre it, also just take the first channel for now
vid_data = squeeze(video(:,:,1,best_startframe:best_endframe));
% we need to shift and rotate the image so that both COR points are 512
hdelta = COR(1)-COR(2);
a = imtranslate(ff, [hdelta/2, 0]);
theta = rad2deg(tan(hdelta/size(video,1)));
a = imrotate(a, theta);
imagesc(a);

vid_data = imrotate(vid_data, theta);
vid_data = imtranslate(vid_data, [hdelta/2, 0]);
%vid_data = imresize(squeeze(video(:,:,1,:)), 0.25);

%% FBP reconstruction
image_width = size(vid_data,2);
num_slices = size(vid_data,1);
% create geometries and projector
proj_geom = astra_create_proj_geom('parallel', 1.0, image_width, angles);
vol_geom = astra_create_vol_geom(image_width,image_width);
recon_id = astra_mex_data2d('create', '-vol', vol_geom);
sinogram_id = astra_mex_data2d('create','-sino', proj_geom);
cfg = astra_struct('FBP');
proj_id = astra_create_projector('linear', proj_geom, vol_geom); % or linear?
cfg.ProjectorId = proj_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ReconstructionDataId = recon_id;
cfg.FilterType = 'Ram-Lak';
cfg.FilterType = 'shepp-logan';

volume = zeros(image_width, image_width, num_slices);
h = waitbar(0, 'Processing FBP');
for i = 1:num_slices-3
    waitbar(i/num_slices, h);
    sino = median(single(vid_data(i:i+2,:,:)),1);% grab a slice through time, median of 3 rows
    sino = squeeze(sino)'; % and transpose it
    astra_mex_data2d('set',sinogram_id,single(sino));
    % run the algorithm
    fbp_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('run', fbp_id);
    volume(:,:,i) = astra_mex_data2d('get', recon_id);
    
end
fprintf('\nDone.\n');
save('reconstructed-mut.mat','volume');
% garbage disposal
astra_mex_data2d('delete', sinogram_id, recon_id);
astra_mex_projector('delete', proj_id);
astra_mex_algorithm('delete', fbp_id);
%%
implay(volume)