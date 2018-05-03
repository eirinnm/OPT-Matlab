function [ best_offset ] = findCOR( sino, angles, leftedge, rightedge, testoffset, focusmethod, plotresults, gpuAvailable )
% tests a bunch of Centers of Rotation
focusfun = @(x) -testCOR(sino, angles, leftedge, rightedge, x, focusmethod, gpuAvailable);
test_offsets = -5+testoffset:5+testoffset;
offset_results = arrayfun(focusfun, test_offsets);
[~, idx] = min(offset_results);
if plotresults
    plot(test_offsets, offset_results);
    title("CoR focus");
end
best_offset = test_offsets(idx);
end