function [ fishmask ] = find_mask( slice_image, thresh)
%FIND_MASK finds the largest object that matches certain shape criteria
SE = strel('square',15); % this is for the dilation step
% bw_slice =  imbinarize(slice_image, thresh);
bw_slice =  imbinarize(slice_image); %auto threshold
cc = bwconncomp(bw_slice, 4);
stats = regionprops(cc, 'Area','Extent','BoundingBox','Eccentricity','Extrema');
% find objects that meet criteria to look like a fish cross section
idx = find([stats.Extent]>0.3 & [stats.Area]>200 & [stats.Eccentricity]<0.95);
% find the biggest object that meets this criteria
fishmask = zeros(size(slice_image));
if ~isempty(idx)
   % we found a suitable object
   [~,i] = max([stats(idx).Area]);
   sliceBB = stats(idx(i)).BoundingBox;
   BW2 = ismember(labelmatrix(cc), idx(i));
   % get the convex hull of this fish slice
   % then dilate it and add it to the mask array
   fishmask = imdilate(bwconvhull(BW2),SE);
end   

end

