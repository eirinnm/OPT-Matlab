videoname = 'Lov-2017-11-06T18_22_43_compressed';
ch1_filename = strcat(videoname,'_ch1.mat');
ch3_filename = strcat(videoname,'_ch3.mat');

% ch1 = load(ch1_filename,'views');
load(ch3_filename,'views');
%%

%%
% bw_views = views>256*level; %make a boolean version, very memory intensive
h=waitbar(0,'Calculating threshold...');
level = graythresh(views);
bw_views = zeros(size(views),'uint8');
numslices = size(bw_views,3);
waitbar(0,h,'Applying threshold and filtering slices...');
for slicenum = 1:numslices
    waitbar(slicenum/numslices,h);
    bw_views(:,:,slicenum) = imbinarize(views(:,:,slicenum),level);
end
close(h);
%%
conncomponents = bwconncomp(bw_views);
%%
numPixels = cellfun(@numel,conncomponents.PixelIdxList);
[biggest,idx] = max(numPixels);
bw_views(:)=0; 
bw_views(conncomponents.PixelIdxList{idx}) = 1;
clear conncomponents
%%
% find the convex contour of the fish in each slice
outlines = zeros(size(views),'uint8');
h=waitbar(0,'Finding outline...');
for slicenum = 1:numslices
    waitbar(slicenum/numslices,h);
    outlines(:,:,slicenum) = bwconvhull(bw_views(:,:,slicenum));
end
close(h);
%%
stats = regionprops(outlines,'BoundingBox');
BB = round(stats.BoundingBox);
%%
% BB = [213   255   142   415   596   883]
cropped_ch3 = views(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1);
%%
clear views
clear bw_views
clear outlines
%%
load(ch1_filename,'views');
%%
cropped_ch1 = views(BB(2):BB(2)+BB(5)-1,BB(1):BB(1)+BB(4)-1,BB(3):BB(3)+BB(6)-1);
%%
outputfilename = strcat(videoname,'_reconstructed_ch1.tiff');
imwrite(cropped_ch1(:,:,1),outputfilename);
for i = 2:size(cropped_ch1,3)
    imwrite(cropped_ch1(:,:,i),outputfilename,'WriteMode', 'append');
end
outputfilename = strcat(videoname,'_reconstructed_ch3.tiff');
imwrite(cropped_ch3(:,:,1),outputfilename);
for i = 2:size(cropped_ch3,3)
    imwrite(cropped_ch3(:,:,i),outputfilename,'WriteMode', 'append');
end
%%

% function [ fishvol ] = find_fish( ch1, ch2 )
% %FIND_FISH Finds the fish in a 3D reconstruction
% 
% 
% 
% end

