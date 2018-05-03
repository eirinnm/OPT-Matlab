function [shifted] = shift_with_nans(img, x)
    shifted = imtranslate(img,[0 x],'nearest');
    shifted(shifted==0)=nan;