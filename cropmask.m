inputname = 'OPT-2018-01-30T16_35_13_recon.mat';
load(inputname);
outline = squeeze(max(mask,[],3));
% mask out the fish
%%
mask = repmat(outline,[1 1 size(mask,3) 3]);
rec(~mask)=0;
clear mask
toppad = max(round((640-size(rec,1))/2),0)
leftpad = max(round((640-size(rec,2))/2),0)
newrec=padarray(rec,[toppad leftpad 0 0]);
clear rec
clear newmask
%%
outputfilename = [inputname(1:end-4),'_cropped.tiff'];
options.overwrite = true;
options.color = true;
options.message = true;
options.compress = 'no';
dat = permute(newrec, [3 1 4 2]);
saveastiff(dat, outputfilename, options);
clear dat