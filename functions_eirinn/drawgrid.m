% draw lines of latitude along the fish
spacing = 50;
meshmatrix = zeros(size(ch3),'uint8');
for slicenum = 50:spacing:size(ch3,3)-10
    I=(squeeze(ch3(:,:,slicenum)));
    BW = imbinarize(I, 0.07);
    BW2 = bwareafilt(BW,1);
    B = bwboundaries(BW2,4,'noholes');
    boundary = B{1};
    x = B{1}(:,2);
    y = B{1}(:,1);
    smoothX = round(sgolayfilt(x, 2, 99));
    smoothY = round(sgolayfilt(y, 2, 99));
    slicenumpoints = ones(size(smoothX)) * slicenum;
    linearInd = sub2ind(size(meshmatrix),smoothY, smoothX, slicenumpoints);
    meshmatrix(linearInd) = 255;
    for nextslice = slicenum:slicenum+5
        meshmatrix(:,:,nextslice) = meshmatrix(:,:,slicenum);
    end
end
for slicenum = 50:spacing:size(ch3,2)-10
    I=(squeeze(ch3(:,slicenum,:)));
    BW = imbinarize(I, 0.07);
    BW2 = bwareafilt(BW,1);
    B = bwboundaries(BW2,4,'noholes');
    boundary = B{1};
    x = B{1}(:,2);
    y = B{1}(:,1);
    smoothX = round(sgolayfilt(x, 2, 99));
    smoothY = round(sgolayfilt(y, 2, 99));
    slicenumpoints = ones(size(smoothX)) * slicenum;
    linearInd = sub2ind(size(meshmatrix),smoothY, slicenumpoints, smoothX);
    meshmatrix(linearInd) = 255;
    for nextslice = slicenum:slicenum+10
        meshmatrix(:,nextslice,:) = meshmatrix(:,slicenum,:);
    end
end

meshmatrix = permute(meshmatrix, [3 1 2]);
options.overwrite=true;
saveastiff(meshmatrix,'meshmatrix.tiff', options);