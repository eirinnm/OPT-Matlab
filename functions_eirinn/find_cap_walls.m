function [ edges, cap_width_scalefactor ] = find_cap_walls(vid_data)
numslices = size(vid_data,1);
capillary_endpoints = find_hough_points(vid_data);
% draw the found capillary walls
found_walls = capillary_endpoints(capillary_endpoints(:,3)==1,1:2);
subplot(2,2,1);
imagesc(vid_data(:,:,1));
title('Frame 1');
for i = 1:size(found_walls)
    vline = line(found_walls(i,:), [1 numslices]);
    set(vline,'color','r');
    set(vline,'LineStyle','--');
    set(vline,'LineWidth',3);
end
% average every row, per frame, to make a blended sinogram
scanzones = squeeze(mean(vid_data,1));
subplot(2,2,2);
imagesc(scanzones);
title("Average sinogram");
% now find the walls in the sinogram
cap_peak_points = find_walls(double(scanzones), capillary_endpoints);

hold on;
plot(cap_peak_points,'color','r','LineStyle',':','LineWidth',3);
hold off
view([-90 -90])
% is the capillary width shifting?
cap_width = diff(cap_peak_points,1,2);
cap_width_smooth = sgolayfilt(cap_width,2,33);
cap_width_mean = mean(cap_width_smooth);
cap_width_scalefactor = (cap_width_smooth-cap_width_mean)/cap_width_mean;
subplot(2,2,4);
plot([cap_width cap_width_smooth]);
title('Capillary width')
variances = [std(cap_peak_points) std(cap_width_smooth)];
disp('StDev of left wall, right wall, and width:');
disp(variances);
edges = round(nanmedian(cap_peak_points));
end