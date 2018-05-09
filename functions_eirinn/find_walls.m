function [ cap_peak_points ] = find_walls( sino, capillary_endpoints )
%FIND_WALLS Looks for likely walls in a sinogram

numframes=size(sino,2);
image_width=size(sino,1);

good_peaks = [];
for framenum = 1:numframes
    smoothslice = sgolayfilt((sino(:,framenum)),3,35);
    capillary_candidates = capillary_endpoints(capillary_endpoints(:,3)==framenum,1);
    % find the peaks in this 1D slice
    [pks, locs] = findpeaks(smoothslice, 'MinPeakProminence', 10,'MinPeakDistance',50, 'SortStr','descend','NPeaks',4,'MaxPeakWidth',70);
    for i=1:size(pks)
        % is this peak nearby a capillary candidate?
        if min(abs(locs(i)-capillary_candidates))<100
            % add this to the list of good peaks
            good_peaks = [good_peaks; pks(i) locs(i) framenum];
        end
    end
end

% Now we have a list of candidate capillary locations. We want to find the
% inner ones. We expect, and require, that the inner ones are visible in every frame.
% But it's possible that the outer ones are also visible in every frame. 
% for each frame, find the 2 candidates that are closest to the image
% centre. These are the inner walls.

% Cluster the peaks:
mean_found_peaks = round(size(good_peaks,1)/numframes);
[good_peaks(:,4), peak_centres] = kmeans(good_peaks(:,2),mean_found_peaks);
% which of these clustered peaks are closest to the centre of the image?
peak_distances = image_width/2 - peak_centres;%
peak_distances(:,2) = peak_distances>0;
peak_distances(:,3) = 1:mean_found_peaks;
left_peaks = peak_distances(peak_distances(:,2)==1,:);
[~, idx] = min(abs(left_peaks(:,1)));
left_peak_idx = left_peaks(idx,3);

right_peaks = peak_distances(peak_distances(:,2)==0,:);
[~, idx] = min(abs(right_peaks(:,1)));
right_peak_idx = right_peaks(idx,3);

% peak_distances(peak_distances>0)
% find2Closest = @(a) {a(mink(abs(a(:,2)-image_width/2),2),:)};
% inner_walls = splitapply(find2Closest,good_peaks,good_peaks(:,3));
% inner_walls = cell2mat(inner_walls);
% 
% cap_wall_left = inner_walls(inner_walls(:,2)<image_width/2,2:3);
% cap_wall_right = inner_walls(inner_walls(:,2)>=image_width/2,2:3);
cap_wall_left = good_peaks(good_peaks(:,4)==left_peak_idx,2:end);
cap_wall_right = good_peaks(good_peaks(:,4)==right_peak_idx,2:end);
% plot([cap_wall_left cap_wall_right]);

% create an empty matrix to hold the actual points per frame
cap_peak_points = zeros(numframes, 2);
cap_peak_points(:) = nan;

cap_peak_points(cap_wall_left(:,2),1) = cap_wall_left(:,1);
cap_peak_points(cap_wall_right(:,2),2) = cap_wall_right(:,1);
% nans are where one or both cap walls aren't found. 
% replace nans with a value based on median cap width.
cap_width = nanmedian(diff(cap_peak_points,1,2));
missing_walls = isnan(cap_peak_points(:,1));
cap_peak_points(missing_walls,1) = cap_peak_points(missing_walls,2)-cap_width;
missing_walls = isnan(cap_peak_points(:,2));
cap_peak_points(missing_walls,2) = cap_peak_points(missing_walls,1)+cap_width;

end



