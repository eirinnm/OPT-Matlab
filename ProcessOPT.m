%%%
% OPT Reconstruction script
% Eirinn Mackay, Wilson lab October 2017
% Inspired by Amin Allalou (Yanik lab)
%%%
% downscalefactor = 2;
use_image_warping = false;
disp('--- OPT Reconstruction script ----------------');
disp('--- Eirinn Mackay, Wilson lab October 2017 ---');
addpath(genpath('OPT-InSitu-Toolbox/optreconstruction'));
addpath(genpath('OPT-InSitu-Toolbox/3rdparty/astra-1.5'));
try
    d = gpuDevice;
    gpuAvailable = d.SupportsDouble;
    disp('GPU available!');
catch
    gpuAvailable = false;
    disp('No Nvidia GPU detected.');
end
cancelled = 0;
[VideoNames,VideoPath,FilterIndex] = uigetfile('*.avi','Please select the AVI file to be reconstructed:','MultiSelect','on');
% if ~iscell(VideoNames) 
%     if VideoNames == 0
%         return
%     end
% end
batchmode=1;
if ischar(VideoNames)
    VideoNames = {VideoNames};
    batchmode=0;
end
for V=1:length(VideoNames) % start a major loop through videos
VideoName = VideoNames{V};
if cancelled == 1
    break
end
%% Read video file and subtract the background.
disp('Reading input video...');
disp(VideoName);
h = waitbar(0.33,'Loading video...');
video = read(VideoReader([VideoPath,VideoName]));
close(h);
% background = imread('background2018-02-23T16_38_32.bmp');
background = imread('background2018-03-09T14_30_27.bmp');
% subtract the background. This will invert the brightness too
video = (background - video);
% How long is the revolution? Find the frames with the least difference
% from the first frame
frame_diffs = squeeze(sum(sum(abs(video(:,:,2,:) - video(:,:,2,1)))));
[~, minidx] = min(frame_diffs(300:end));
minidx = minidx + 300;
fprintf("Estimated rotation length: %d frames\n",minidx);
% clip the video
video = video(:,:,:,1:minidx);
numslices = size(video,1);
image_width = size(video,2);
numchannels = size(video,3);
numframes = size(video,4);
% create a figure to show progress
f = figure;

%% Find capillary walls
[~, cap_width_scalefactor, cap_peak_points] = find_cap_walls(squeeze(video(:,:,2,:)));
% which frame has the capillary closest to the centre of the image?
cap_mid_distance_to_centre = mean(cap_peak_points,2)-image_width/2;
[~, closest_frame] = min(abs(cap_mid_distance_to_centre));
initial_shift_by = cap_mid_distance_to_centre(closest_frame);

%% Video stabilization
subplot(2,2,1);
imagesc(video(:,:,2,1));
title('Frame 1');
scanzones = squeeze(mean(video(:,:,2,:),1));
subplot(2,2,2);
imagesc(scanzones');
title("Average sinogram");

% for each consecutive frame, find the lateral movement that minimizes the
% difference (least squares)
shift_by = stabilize_capillary(scanzones, initial_shift_by, closest_frame);
shift_by = cumsum(shift_by);
%new: centre the shift amount to the position from the frame closest to the
%middle of the image
shift_by = shift_by-shift_by(closest_frame);%+initial_shift_by;
scanzones_fixed = scanzones;
for framenum = 1:numframes
    scanzones_fixed(:,framenum) = imtranslate(scanzones(:,framenum),[0 shift_by(framenum)],'linear');
end
%where do we think the capillary walls are now?
adjusted_cap_walls = cap_peak_points + shift_by; 
approx_cap_wall = round(mean(adjusted_cap_walls));
% Note! Not adjusted for capillary width variance
subplot(2,2,3);
imagesc(scanzones_fixed);
title('Stabilized sinogram');
hold on;
plot(adjusted_cap_walls,'color','r','LineStyle',':','LineWidth',3);
hold off
view([-90 -90])

%% Adjust the entire video
% Resize the video horizontally so that the capillary is the same width in all frames.
% Also shift the video horizontally to stabilise it.
h=waitbar(0,'Stabilizing and centering video...');
image_width = size(video,2);
image_height = size(video,1);
% video2 = video;
for framenum = 1:numframes
    I = video(:,:,:,framenum);
    %find the desired width of this frame
    newwidth = round(image_width*(cap_width_scalefactor(framenum)*-1+1));
    %where is the new centre of that image in relation to the old centre?
    offset_to_centre = round((image_width - newwidth)/2);
    % We also need to move the video by the amount given in shift_by
    offset_to_centre = offset_to_centre + shift_by(framenum);
    % Actual transformations:
    % resize the video (this will resize from the left side)
    Ir = imresize(I,[image_height newwidth]);
    % now shift it so it's in the centre
    Ir = imtranslate(Ir,[offset_to_centre 0]);
    % copy it back into the video matrix
    newwidth_int = floor(newwidth);
    if newwidth_int<image_width
        Ir = Ir(1:image_height,1:newwidth_int,:);
        video(:,1:newwidth_int,:,framenum) = Ir;
    else
        Ir = Ir(1:image_height,1:image_width,:);
        video(:,:,:,framenum) = Ir;
    end
    waitbar(framenum/numframes,h);
end
close(h);

%crop the video to the stabilized region
%note: this will crash if the cap walls are closer than 50px from the image
%edge
crop_left = approx_cap_wall(1)-50;
crop_right = approx_cap_wall(2)+50;
video = video(:,crop_left:crop_right,:,:);
% Find walls in the stabilized video
[edges, cap_width_scalefactor, cap_peak_points] = find_cap_walls(squeeze(video(:,:,2,:)));
% these edges will be used to crop and centre the sinogram prior to
% reconstruction

%% Reconstruct a slice
ch = 2; %channel 2 is the green channel which has a mix of ISH and SYTOX
slicenum = 700; % this is the slice with the highest intensity
subplot(2,2,1);
imagesc(video(:,:,ch,1));
hline = refline(0,slicenum);
    set(hline,'color','r');
    set(hline,'LineStyle','--');
title(sprintf('Test slice (ch%d)', ch));
line([edges(1) edges(1)],[1 numslices], 'color','r');
line([edges(2) edges(2)],[1 numslices], 'color','r');

angles = linspace(0,2*pi, numframes);

sino = squeeze(video(slicenum,:,ch,:));
sino(~imbinarize(sino))=0;
% imagesc(sino)
contract_borders_by = -20;
leftedge = edges(1)+contract_borders_by;
rightedge = edges(2)-contract_borders_by;
plotresults = true;
subplot(2,2,2);
offset_search_region = 0;
extra_offset = findCOR(sino, angles, leftedge,rightedge,offset_search_region,'GDER',plotresults,gpuAvailable);
fprintf('Optimised COR offset is %0.2f pixels\n',extra_offset);

extra_offset=extra_offset+0;
% display this reconstruction
contract_borders_by = -20;
leftedge = edges(1)+contract_borders_by+extra_offset;
rightedge = edges(2)-contract_borders_by+extra_offset;
subplot(2,2,3);
fixed_sino = sino(leftedge:rightedge,:);
imagesc(fixed_sino');
title('Sinogram for reconstruction');
% adjust the edges again to remove the capillary. This will get used in the
% final reconstruction.
contract_borders_by = +40;
leftedge = edges(1)+contract_borders_by+extra_offset;
rightedge = edges(2)-contract_borders_by+extra_offset;

subplot(2,2,4);
downsamplefactor = 1;
small_sino = imresize(fixed_sino, 1/downsamplefactor);
small_angles = decimate(angles,downsamplefactor);
rec_slice = reconstruct_from_sino(small_sino', small_angles,gpuAvailable);
rec_slice(rec_slice<0)=0;
imagesc(rec_slice);
title('Reconstructed slice');
daspect([1 1 1]);

%% Reconstruct all slices if it looks good
clear Movie;
framenum=1;
chans = [2 1 3];
image_width = size(fixed_sino,1);
if batchmode==0
    button = questdlg('Proceed with full reconstruction?','Hey','Yes','No','Yes');
    cancelled=0;
    if isequal(button,'No')
        cancelled = 1;
    end
end
if cancelled==0 %okay let's do this
    delete(f);
    downsamplefactor = 1;
    display_subplots = false;
    vol_width = ceil(image_width/downsamplefactor);
    vol_height = ceil(numslices/downsamplefactor);
    reconTimer = tic;
    vol_width = rightedge-leftedge+1;
    views = zeros(vol_width, vol_width, 3, vol_height,'uint8'); % this will hold our data
    fishmask = zeros(vol_width, vol_width, vol_height,'logical');
    for ch = chans
        if cancelled == 1
            break
        end
        fprintf('Reconstructing channel %d\n', ch);
        max_projection = zeros(vol_width, vol_height, 'uint8');
        first_frame = video(:,:,ch,1); %purely for subplot illustration
        h = waitbar(0,'Processing...','CreateCancelBtn','cancelled = 1; delete(h)');
        f = figure('CloseRequestFcn','cancelled = 1; delete(f);');
        slice_time_elapsed = 0.4; % start with a reasonable number for time estimation
        for slicenum = 1:downsamplefactor:numslices
            % check for cancel button
            if cancelled == 1
                break
            end
            tic
            seconds_remaining = slice_time_elapsed*(numslices - slicenum)/downsamplefactor;
            hms = fix(mod(seconds_remaining, [0, 3600, 60]) ./ [3600, 60, 1]);
%             fprintf('Time remaining: %02d:%02d:%02d\n', hms);
%             waitbar(slicenum/numslices, h, sprintf('Processing slice %d',slicenum));
            waitbar(slicenum/numslices, h, sprintf('Time remaining: %02d:%02d:%02d', hms));
            sino = (squeeze(video(slicenum,:,ch,:)));
            fixed_sino = sino(leftedge:rightedge,:);
            small_sino = imresize(fixed_sino, 1/downsamplefactor);
            small_angles = decimate(angles,downsamplefactor);
            % actually do the reconstruction
            rec_slice = reconstruct_from_sino(small_sino', small_angles,gpuAvailable);
            rec_slice(rec_slice<0)=0;
            rec_slice = uint8(round(rec_slice*64));
            % filter it
            rec_slice = medfilt2(rec_slice);
            % add the reconstruction to the volume
            downsampled_slicenum = uint16((slicenum-1)/downsamplefactor+1);
            views(:,:,ch,downsampled_slicenum) = rec_slice;
            % if this is channel 3, find the fish outline
            if ch == 2
                this_fishmask = find_mask(rec_slice, 0.07);
                fishmask(:,:,downsampled_slicenum) = this_fishmask;
            end
            slice_projection = max(rec_slice,[],1);
            max_projection(:,downsampled_slicenum) = slice_projection;
            
            if display_subplots
                subplot(2,2,1);
                imagesc(first_frame);
                hline = refline(0,slicenum);
                set(hline,'color','r');
                set(hline,'LineStyle','--');
                title(sprintf('Ch %d slice %d',ch,slicenum));
                subplot(2,2,2);
                imagesc(small_sino')
                title('Sinogram');
                subplot(2,2,3)
                imagesc(rec_slice)
                daspect([1 1 1]);
                title('Reconstruction (top view)')
%                 if ch == 3
%                     % draw the bounding box
%                     max_mask = squeeze(max(fishmask,[],3));
%                     whole_region = regionprops(max_mask);
%                     if ~isempty(whole_region)
%                         wholeBB = whole_region.BoundingBox;
%                         hold on
%                         rectangle('Position',wholeBB, 'EdgeColor','green','LineStyle','--');
%                         hold off
%                     end
%                 end
                subplot(2,2,4)
                imagesc(max_projection');
                daspect([1 1 1]);
                title('Max intensity projection')
                Movie(framenum) = getframe(gcf);
                framenum=framenum+1;
            else
                subplot(1,2,1);
                imagesc(rec_slice);
                daspect([1 1 1]);
                title(sprintf('Processing slice %d of %d',slicenum,numslices));
                subplot(1,2,2);

                imagesc(max_projection');
                daspect([1 1 1]);
                title('Projection')
            end
                    
            slice_time_elapsed = toc;
        end

        delete(f);
        delete(h);
    end
    % all channels complete
    hms = fix(mod(toc(reconTimer), [0, 3600, 60]) ./ [3600, 60, 1]);
    disp(sprintf("Elapsed time: %02d:%02d:%02d", hms));
    if cancelled == true
        disp('Processing cancelled.')
    else
        disp('Saving cropped volume...');
        h=waitbar(1,'Saving cropped volume to file...');
        
        %Crop the volume to the fish
        % Let's just take a max through the
        % stack and find the biggest object there.
        max_mask_topdown = squeeze(max(fishmask,[],3));
        max_mask_props = regionprops(max_mask_topdown);
        [~,mask_biggest_idx] = max([max_mask_props.Area]);
        BB = round(max_mask_props(mask_biggest_idx).BoundingBox);
        rec = views(BB(2):BB(2)+BB(4)-1, BB(1):BB(1)+BB(3)-1, :, :);
        mask = fishmask(BB(2):BB(2)+BB(4)-1, BB(1):BB(1)+BB(3)-1, :);
        save([VideoPath,VideoName(1:end-4),'_recon.mat'],'rec', 'mask','-v7.3');
        % Save the output to a TIFF stack
        outputfilename = [VideoPath,VideoName(1:end-4),'_recon.tiff'];
        options.overwrite = true;
        options.color = true;
        options.message = false;
        saveastiff(rec, outputfilename, options);
        disp('Done.');
        delete(h);
        clear dat
    end
end

end %% end the major loop through videos

%% Attempt a 3d scatter plot render
% [xx,yy,zz] = meshgrid(1:134,1:119,1:256);
% n = nonzeros(Ds(Ds>5));
% n = find(D>5);
% [xx,yy,zz] = ind2sub(size(D),n);
% scatter3(xx,yy,zz,1,'Marker','.','MarkerEdgeAlpha',0.01);
% axis vis3d
% daspect([1,1,1]);
%% Prepare a downsampled volume mask
% meanrec = mean(rec,4);
% % B = imresize3(meanrec,0.25);
% B = imresize3(rec(:,:,:,2),0.25);
% BW = B>(10);
% CC = bwconncomp(BW);
% numOfPixels = cellfun(@numel,CC.PixelIdxList);
% [unused,indexOfMax] = max(numOfPixels);
% biggest = zeros(size(B),'logical');
% biggest(CC.PixelIdxList{indexOfMax}) = 1;
% %%
% BW = B>3;
% imagesc(squeeze(BW(80,:,:)))