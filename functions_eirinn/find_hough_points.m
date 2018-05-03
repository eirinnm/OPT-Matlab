function [ capillary_endpoints ] = find_hough_points( vid_data )
%FIND_HOUGH_HINTS uses Hough transform to find vertical lines
% Use hough transform to find lines
resizefactor = 4;
numframes = size(vid_data,3);
capillary_endpoints = [];
h = waitbar(0,'Finding capillary walls with Hough transform');
for framenum = 1:numframes
    waitbar(framenum/numframes);
    a = vid_data(:,:,framenum);
    a = imresize(a, 1/resizefactor);
    I = edge(a, 'Sobel');
    % Hough transform
    [H,T,R]= hough(double(I),'Theta',-1:0.5:1); % Find vertical lines by limiting theta to a small range
    P  = houghpeaks(H,4);%
    lines = houghlines(I,T,R,P);
    % extend all the lines to the top and bottom
    for l = 1:size(lines,2)
        y=[lines(l).point1(2) lines(l).point2(2)];
        x=[lines(l).point1(1) lines(l).point2(1)];
        capillary_endpoints = [capillary_endpoints; interp1(y, x, [1 size(a,1)], 'linear', 'extrap')*resizefactor framenum];
    end
end
close(h);

end

