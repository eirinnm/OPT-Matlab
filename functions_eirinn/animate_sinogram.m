%% Animate making the sinogram, for presentations
% 
slicenum=600;
sino=zeros(numframes,image_width);
for framenum = 1:numframes
    subplot(2,1,1);
    thisframe = video(:,:,ch,framenum);
    imagesc(thisframe,[1 127]);
    hline = refline(0,slicenum);
    set(hline,'color','r');
    set(hline,'LineStyle','--');
    daspect([1 1 1]);
    sino(framenum,:) = thisframe(slicenum,:);
    subplot(2,1,2);
    imagesc(sino,[1 127]);
    daspect([2.5 1 1])
%     pause(0.05);
    Movie(framenum) = getframe(gcf);
end
myVideo = VideoWriter([VideoPath,VideoName(1:end-4),'_sinogram.mp4'],'MPEG-4');
open(myVideo);
writeVideo(myVideo,Movie);
close(myVideo);