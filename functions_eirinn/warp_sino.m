function [ shifted_sino ] = warp_sino( sino, offset )
%WARP_SINO Offsets each pixel column in a sinogram
% Sino: 2D image. Offset: array of vertical offset amounts
dfield = repmat(-offset,1,size(sino,1))';
dfield(:,:,2) = dfield(:,:,1);
dfield(:,:,1) = 0;
shifted_sino = imwarp(sino, dfield);
% imagesc(shifted_sino);
end

