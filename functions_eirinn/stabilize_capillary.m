function [ shift_by ] = stabilize_capillary(scanzones)
%%% Steps through lines of the "scanzone" (sinogram of all slices) and
%%% calculates horizontal pixel offsets that minimizes the difference. 

    function [fval] = shift_compare(x)
        translated = imtranslate(scanline_next,[0 x],'linear');
        translated(translated==0) = nan;
        fval = nansum((scanline - translated).^2);
    end
    numframes = size(scanzones,2);
    shift_by = zeros(numframes,1);
%     shift_by = NaN(numframes,1);
%     options = optimset('TolX',0.1,'TolFun',0.1,'Display','iter');
    h=waitbar(0,'Measuring capillary movement...');
    scansize=4;
    for framenum = 1:scansize:numframes-scansize
%         framenum = 361;
        scanline = scanzones(:,framenum);
%         scanline = scanzones(:,150);
        scanline_next = scanzones(:,framenum+scansize);
        test_shifts = -10:10;
        shift_results = arrayfun(@shift_compare, test_shifts);
        [~, idx] = min(shift_results);
        best_shift = test_shifts(idx);
%         shift_by(framenum+scansize) = best_shift;
        shift_by (framenum+1:framenum+scansize) = best_shift/scansize;
        waitbar(framenum/numframes,h);
    end
    close(h);
    
%%% Pixel offsets are in whole-pixel units which means rounding errors
%%% accumulate over time. Run a second pass on every 4th row
end



    