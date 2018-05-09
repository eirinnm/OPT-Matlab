function [ shift_by ] = stabilize_capillary(scanzones, initial_shift_by, closest_frame)
%%% Steps through lines of the "scanzone" (sinogram of all slices) and
%%% calculates horizontal pixel offsets that minimizes the difference. 
%%% New: start at the frame where the capillary is closest to the middle.
    function [fval] = shift_compare(x)
        translated = imtranslate(scanline_next,[0 x],'linear');
        translated(translated==0) = nan;
        fval = nansum((scanline - translated).^2);
    end
    numframes = size(scanzones,2);
    shift_by = zeros(numframes,1);
    h=waitbar(0,'Measuring capillary movement...');
    scansize=4; %we actually look at every 4th frame
    for framenum = 1:scansize:numframes-scansize
        scanline = scanzones(:,framenum);
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
    
end



    