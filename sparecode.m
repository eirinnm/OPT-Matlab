Test this on a sinogram
slicenum = 900;
sino = double(squeeze(vid_data(slicenum,:,:)) -0);
% imagesc(sino);
% title('Uncorrect sinogram');
% %get the capillary middle
% hold on;
% raw_middle = cap_mid(slicenum,:);
% plot(raw_middle,'-', 'LineWidth',3,'Color','green');
% fitted_middle = sinogram_offsets(slicenum,:);
% corrected_middle = (fitted_middle+raw_middle);
% plot(corrected_middle,'--', 'LineWidth',3,'Color','magenta');
% hold off
s = warp_sino(sino, offsetgrid(slicenum,:));
imagesc(s);


% cap_peak_points = zeros(numframes, 2);
% tic
% for framenum = 1:numframes
%     smoothslice = sgolayfilt((sino(:,framenum)),3,55);
%     % find the highest peak within 200 pixels of the left side
%     [pks, locs] = findpeaks(smoothslice(1:200), ...
%         'MinPeakProminence', 20,'MinPeakDistance',50, 'SortStr','descend','NPeaks',1);
%     if isempty(locs)
%         cap_peak_points(framenum,1) = nan;
%     else
%         cap_peak_points(framenum,1) = locs;
%     endplot(top_midpoints, '--', 'LineWidth',3,'Color','green');
%     % find the highest peak within 200 pixels of the right side
%     [pks, locs] = findpeaks(smoothslice(imagewidth-200:imagewidth), ...
%         'MinPeakProminence', 20,'MinPeakDistance',50, 'SortStr','descend','NPeaks',1);
%     if isempty(locs)
%         cap_peak_points(framenum,2) = nan;
%     else
%         cap_peak_points(framenum,2) = locs+imagewidth-200;
%     end
% end
% toc
% wall_finding_rate = (numframes-sum(isnan(cap_peak_points)))/numframes*100
% % find the distance between the walls when the walls are found
% cap_peak_points = medfilt1(cap_peak_points);
% cap_width_points = (diff(cap_peak_points,1,2));
% cap_width = nanmedian(cap_width_points)
% cap_width_stdev = nanstd(cap_width_points)
% [m r] = max(wall_finding_rate);
% if m<90
%     error('Walls not found in enough frames');
% end
% %otherwise, plot the centre point. Either add it or subtract it from the known wall
% if r==1
%     midpoints=cap_peak_points(:,r)+cap_width/2;
% else
%     midpoints=cap_peak_points(:,r)-cap_width/2;
% end
% hold on
% plot(cap_peak_points, '--', 'LineWidth',3,'Color','red');
% plot(midpoints, '--', 'LineWidth',3,'Color','green');
% hold off

% imagesc(sino);
% sinfit = sin_fit(midpoints); %fit to a sine curve, returns amplitude, period, phase, offset
% phase = sinfit(2)/(2*sinfit(3)); 
% fitted = sin(-deg2rad(([1:size(midpoints,1)])*360/sinfit(2) + phase))*sinfit(1);
% fitted = fitted + size(sino/2); %shift the sine wave to be in the centre
% fitted = fitted';
% offset = fitted-midpoints; %how much do we need to warp the sinogram?
% offset = fillmissing(offset, 'linear');
% hold on
% plot(fitted, '-', 'LineWidth',3,'Color','green')
% hold off

Find capillary walls in the first frame
% a = vid_data(:,:,1);
% numslices = size(vid_data,1);
% imagewidth = size(vid_data,2);
% numframes = size(vid_data,3);
% imagesc(a);
% daspect([1 1 1])
% title('Finding capillary walls');
% hold on;
% 
% I = edge(imgaussfilt(a,3), 'Sobel', []);
% 
% % Hough transform
% [H,T,R]= hough(double(I)); % Find the 2 strongest peaks
% P  = houghpeaks(H,3);%
% 
% lines = houghlines(I,T,R,P);
% % get the X coordinate for these lines at every row
% capillary_starts = [];
% for l = 1:size(lines,2)
%     y=[lines(l).point1(2) lines(l).point2(2)];
%     x=[lines(l).point1(1) lines(l).point2(1)];
%     capillary_starts = [capillary_starts; interp1(y, x, [1:numslices], 'linear', 'extrap')];
%     scatter(x,y);
% end
% capillary_starts = sort(capillary_starts', 2);
% % sometimes we end up with multiple candidates for each of the two capillary walls
% % we need to merge nearby points
% merged_cap_starts=zeros(size(capillary_starts,1),2);
% if size(capillary_starts,2)>2
%     clusters = kmeans(capillary_starts(1,:)',2);
%     merged_cap_starts(:,1) = mean(capillary_starts(:,clusters==1),2);
%     merged_cap_starts(:,2) = mean(capillary_starts(:,clusters==2),2);
% end
% capillary_starts = sort(merged_cap_starts,2);
% plot(capillary_starts, [1:numslices],'--', 'LineWidth',3,'Color','red');
% hold off;
% 
% %%  To do: repeat this every 100 frames to get some snapshots of the ideal capillary position

Find capillary walls by looking for peaks in 1D slices
% a = vid_data(:,:,250);
% numslices = size(vid_data,1);
% image_width = size(vid_data,2);
% numframes = size(vid_data,3);
% imagesc(a);
% daspect([1 1 1])
% title('Finding capillary walls');
% % hold on;
% thisframe = imgaussfilt(a,3);
% [toppeaks, toplocs] = findpeaks(double(thisframe(1,:)), 'MinPeakProminence', 20,'MinPeakDistance',50);
% [bottompeaks, bottomlocs] = findpeaks(double(thisframe(end,:)), 'MinPeakProminence', 20,'MinPeakDistance',50);
% % is there a peak on the right? find the highest peak within 300 pixels of the right edge
% % rightside = toppeaks(toplocs > image_width-2000)
% % [m, i] = max(toppeaks(toplocs > image_width-2000))
% % topright = toplocs(toplocs > image_width-2000)(i)
% % [m, i] = max(toppeaks(toplocs > image_width-200));
% % rowpeaks = [];
% % rowlocs = [];
% % allpeaks = [];
% % alllocs = [];
% % for slicenum = 1:300:numslices
% %     thisslice = double(thisframe(slicenum,:));
% %     [pks, locs] = findpeaks(thisslice, 'MinPeakProminence', 20,'MinPeakDistance',50);
% %     rowpeaks=[rowpeaks, pks];
% %     rowlocs=[rowlocs, locs];
% %     allpeaks=[allpeaks, pks];
% %     alllocs=[alllocs, locs];
% %     
% %     %what's the highest peak on the right side of the image?
% % %     disp(size(locs));
% % end
% 
% 
% hold off;

Find capillary walls by looking for peaks at the top
% thisframe = vid_data(:,:,1);
% thisframe = imgaussfilt(thisframe,3);
% imagesc(thisframe);
% daspect([1 1 1])
% slicenum = 1;
% thisslice = double(thisframe(slicenum,:));
% findpeaks(thisslice, 'MinPeakProminence', 50,'MinPeakDistance',50);
% [pks, locs] = findpeaks(thisslice, 'MinPeakProminence', 50,'MinPeakDistance',50);
% % have we found one on the left and right?
% capillary_start = capillary_starts(1,:);
% distances = abs(locs - capillary_start);
% [d, rank] = min(distances)
% if max(d)>50
%     error('No peaks found near expected capillary wall');
% end
% cap1 = locs(rank(1))
% cap2 = locs(rank(2))


% search for the walls across a lot of sinograms
% cap_walls = zeros(numslices,numframes,2);
% search_range = 15;
% slice_errors = zeros(numslices,2);
% tic;
% h = waitbar(0,'Finding capillary walls');
% for slicenum = 1:50:numslices
%     sino = double(squeeze(vid_data(slicenum,:,:)) -0);
%     waitbar(slicenum/numslices,h);
%     last_cap1 = capillary_starts(slicenum,1);
%     last_cap2 = capillary_starts(slicenum,2);
% %     cap_walls(slicenum,1,1) = last_cap1;
% %     cap_walls(slicenum,1,2) = last_cap2;
%     for x = 1:2:numframes %check every 2nd column in the sinogram
%         smoothslice = sgolayfilt((sino(:,x)),3,55);
%         [pks, locs] = findpeaks(smoothslice, 'MinPeakProminence', 10,'MinPeakDistance',50);
%         % find nearest cap1
%         [d, rank] = min(abs(locs-last_cap1));
%         if d < search_range
%             last_cap1 = locs(rank);
%             cap_walls(slicenum,x,1) = last_cap1;
%         else
%             cap_walls(slicenum,x,1) = nan;
%     %         break
%         end
%         % find nearest cap2
%         [d, rank] = min(abs(locs-last_cap2));
%         if d < search_range
%             last_cap2 = locs(rank);
%             cap_walls(slicenum,x,2) = last_cap2;
%         else
%             cap_walls(slicenum,x,2) = nan;
%     %         break
%         end
%     end
%     slice_errors(slicenum,1) = max(abs(diff(cap_walls(slicenum,:,1))));
% end
% close(h);
% % cap_mid = mean(cap_walls,3);
% % cap_mid(cap_mid==0)=nan;
% cap_walls(cap_walls==0)=nan;
% cap_walls = fillmissing(cap_walls, 'linear', 1); % fill the skipped slices in the stack
% cap_walls = fillmissing(cap_walls, 'linear', 2); % fill the gaps in the sinogram
% cap_mid = mean(cap_walls,3);
% toc

Fit the midpoint to ideal sine curves
% sinogram_offsets = zeros(numslices, numframes);
% real_sinogram_offsets = [];
% for slicenum = 1:50:numslices
%     midpoints = cap_mid(slicenum,:)'; %this is the midpoint of the capillary sinogram
%     y = midpoints;
%         x = [1:size(y)]';
%         yu = max(y);
%         yl = min(y);
%         yr = (yu-yl);                               % Range of ‘y’
%         yz = y-yu+(yr/2);  % y centred around the origin
%         z = yz .* circshift(yz,[0 1]);
%         zx = x(z <= 10);     % Find zero-crossings
%         zdiff = diff(zx);
%         per = 2*mean(zdiff(zdiff>1));                     % Estimate period
%         ym = mean(y);                               % Estimate offset
%         
%         fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);    % Function to fit
%         fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
%         sinfit = fminsearch(fcn, [yr;  per;  -1;  ym]);                    % Minimise Least-Squares
% 
%     phase = sinfit(2)*sinfit(3); % not sure about this
% %     phase = 180;
%     fitted = sin(-deg2rad(([1:size(midpoints,1)])*360/sinfit(2)-phase+180))*sinfit(1)+512;
%     fitted = fitted';
%     offsets = fitted-midpoints;
%     offsets = fillmissing(offsets, 'linear')';
%     real_sinogram_offsets = [real_sinogram_offsets; offsets];
% end
% surf(real_sinogram_offsets, 'edgecolor','none');
% view(2);
% good_slices = [1 5 8 17 20 22];
% sinogram_offsets(good_slices*50-49,:) = real_sinogram_offsets(good_slices,:);
% sinogram_offsets(sinogram_offsets==0)=nan;
% sinogram_offsets = fillmissing(sinogram_offsets,'linear',1);
% clipping_walls = zeros(numslices, numframes,2);
% clipping_walls(good_slices*50-49,:) = cap_walls(good_slices,:);
% clipping_walls(clipping_walls==0)=nan;
% clipping_walls = fillmissing(clipping_walls,'linear',1);
% surf(sinogram_offsets, 'edgecolor','none');
% view(2);

What's the goodness of fit?
look vertically as it's nearly linear that way
% cap_walls_fixed = zeros(size(cap_walls));
% for framenum = 1:2:numframes
%     y=cap_walls(:,framenum,1);
%     y=y(~isnan(y));
%     X=1:numslices;
%     X=X(~isnan(y));
%     fitmode = 'poly1';
%     f1 = fit(X',y,fitmode);
%     %find outliers
%     fdata = feval(f1,X);
%     I = abs(fdata - y) > std(y); %% any outliers?
%     if sum(I)>0
%         outliers = excludedata(X',y,'indices',I);
%         y(outliers)=nan;
%         f2 = fit(X',y,fitmode,'Exclude',outliers);
%         fdata2 = feval(f2,X);
%         I = abs(fdata2 - y) > 0.2*nanstd(y);
%         y(I) = nan;
%         y = fillmissing(y,'linear');
%     end
%     cap_walls_fixed(X,framenum,1) = y;
% end

% cap_walls = fillmissing(cap_walls, 'linear', 1); % fill the skipped slices in the stack
% cap_walls = fillmissing(cap_walls, 'linear', 2); % fill the gaps in the sinogram
% surf(cap_walls(:,:,1), 'edgecolor','none');
% view(2);
% title('Before interpolation');
% cap_walls_fixed(cap_walls_fixed==0)=nan;
% cap_walls_fixed = fillmissing(cap_walls_fixed, 'linear', 1); % fill the skipped slices in the stack
% cap_walls_fixed = fillmissing(cap_walls_fixed, 'linear', 2); % fill the gaps in the sinogram
% figure;
% surface(cap_walls_fixed, 'edgecolor','none');
% title('After interpolation');

% 
% findpeaks(smoothslice, 'MinPeakProminence', 10,'MinPeakDistance',50, 'MaxPeakWidth',40)   ;
% allpeaks = [];
% for x = 1:size(sino,2)
%     smoothslice = sgolayfilt(double(sino(:,x)),5,55);
%     [pks, locs] = findpeaks(smoothslice, 'MinPeakProminence', 10,'MinPeakDistance',50)   ;
%     allpeaks = [allpeaks; {pks, locs}];
% end
% p = vertcat(allpeaks{:});
% % allpeaks = sort(allpeaks,2);
% % title('Capillary wall peak-finding');
% % xlabel('X coordinate');
% % ylabel('Pixel intensity');




% function peaks = get_capillary_walls(sinogram)
%     peaks = [];
%     for t = 1:size(sinogram,2)
%         smoothslice = sgolayfilt(double(sinogram(:,t)),5,55);
%         [pks,loc] = findpeaks(smoothslice, 'MinPeakProminence', 10,'MinPeakDistance',300, 'SortStr','descend');
%         % get the first two peaks by height and average their position
%         if size(loc,1)>=2
%             peaks = [peaks; loc(1:2)'];
%         else
%             peaks = [peaks; [nan, nan]];
%         end
%     end
%     peaks = sort(peaks, 2);
% end

% function midpoints = get_capillary_midpoint(sino)
%     numframes=size(sino,2);
%     imagewidth=size(sino,1);
%     cap_peak_points = zeros(numframes, 2);
%     tic
%     for framenum = 1:numframes
%         smoothslice = sgolayfilt((sino(:,framenum)),3,55);
%         % find the highest peak within 200 pixels of the left side
%         [pks, locs] = findpeaks(smoothslice(1:200), ...
%             'MinPeakProminence', 20,'MinPeakDistance',50, 'SortStr','descend','NPeaks',1);
%         if isempty(locs)
%             cap_peak_points(framenum,1) = nan;
%         else
%             cap_peak_points(framenum,1) = locs;
%         end
%         % find the highest peak within 200 pixels of the right side
%         [pks, locs] = findpeaks(smoothslice(imagewidth-200:imagewidth), ...
%             'MinPeakProminence', 20,'MinPeakDistance',50, 'SortStr','descend','NPeaks',1);
%         if isempty(locs)
%             cap_peak_points(framenum,2) = nan;
%         else
%             cap_peak_points(framenum,2) = locs+imagewidth-200;
%         end
%     end
%     toc
%     wall_finding_rate = (numframes-sum(isnan(cap_peak_points)))/numframes*100
%     % find the distance between the walls when the walls are found
%     cap_peak_points = medfilt1(cap_peak_points);
%     cap_width_points = (diff(cap_peak_points,1,2));
%     cap_width = nanmedian(cap_width_points)
%     cap_width_stdev = nanstd(cap_width_points)
%     [m r] = max(wall_finding_rate);
%     if m<90
%         error('Walls not found in enough frames');
%     end
%     %otherwise, plot the centre point. Either add it or subtract it from the known wall
%     if r==1
%         midpoints=cap_peak_points(:,r)+cap_width/2;
%     else
%         midpoints=cap_peak_points(:,r)-cap_width/2;
%     end
% end
% 
% 
% 
% function [offset, midpoints, fitted] = get_sin_offset(sino)
%     midpoints = get_capillary_midpoint(sino);
%     sinfit = sin_fit(midpoints); %fit to a sine curve, returns amplitude, period, phase, offset 
%     phase=sinfit(3);
%     fitted = sin(deg2rad(phase + ([1:size(midpoints,1)])*360/sinfit(2)))*sinfit(1);
% %     fitted = sin(-deg2rad(([1:size(midpoints,1)])*360/sinfit(2) + phase))*sinfit(1);
% %     fitted = -sin(-deg2rad(phase - ([1:size(midpoints,1)])*360/sinfit(2)))*sinfit(1);
%     fitted = fitted + size(sino,1)/2; %shift the sine wave to be in the centre
%     fitted = fitted';
%     offset = fitted-midpoints; %how much do we need to warp the sinogram?
%     offset = fillmissing(offset, 'linear');
% end
% 
% function S = warp_sino(sino, offset)
%     dfield = ones(size(sino));
%     dfield = dfield(:,1) .* -offset;
%     dfield(:,:,2) = dfield(:,:,1);
%     dfield(:,:,1) = 0;
%     shifted_sino = imwarp(sino, dfield);
%     S = shifted_sino;
% end