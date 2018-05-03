%% Animate aligning the CoR
slicenum=600;
midpoint = mean(edges);
sino = squeeze(video(slicenum,:,ch,:));
targetoffset = image_width/2 - midpoint;
clear Movie
framenum=1;
% offset = -48;
for offset = 0:-1:-48
    fixed_sino = imtranslate(sino, [0 offset]);
    subplot(1,2,1);
    imagesc(fixed_sino);
    hold on;
    midpoint = mean(edges) + offset;
    vline = line([1 numframes], [midpoint midpoint]);
    set(vline,'color','r');
    set(vline,'LineStyle','--');
    set(vline,'LineWidth',1);
    hold off
    view([-90 -90])
    daspect([1 3 1]);
    rec_slice = reconstruct_from_sino(fixed_sino', angles,gpuAvailable);
    rec_slice(rec_slice<0)=0;
    subplot(1,2,2);
    imagesc(rec_slice);
    daspect([1 1 1]);
    Movie(framenum) = getframe(gcf);
    framenum=framenum+1;
end
myVideo = VideoWriter([VideoPath,VideoName(1:end-4),'_alignCOR.mp4'],'MPEG-4');
myVideo.FrameRate=10;
open(myVideo);
writeVideo(myVideo,Movie);
close(myVideo);