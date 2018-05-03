function [ image_focus ] = testCOR( sino, angles, leftedge, rightedge, testoffset, focusmethod, gpuAvailable )
%TESTCOR reconstruct and find the focus value for a given offset
downsamplefactor = 2;
% fixed_sino = warp_sino(sino, knownoffset+testoffset);
% fixed_sino = imtranslate(sino,[0 testoffset]);
testoffset = round(testoffset);
fixed_sino = sino(leftedge+testoffset:rightedge+testoffset,:);
small_sino = imresize(fixed_sino, 1/downsamplefactor);
small_angles = decimate(angles,downsamplefactor);
view = reconstruct_from_sino(small_sino', small_angles,gpuAvailable);
view(view<0)=0;
image_focus = fmeasure(view,focusmethod);

end

