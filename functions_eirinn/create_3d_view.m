%% Hard mask parts of the volume that aren't in the top-down mask
max_mask_topdown = squeeze(max(mask,[],3));
rec_masked = rec;
numslices = size(rec_masked,4);
for slicenum=1:numslices
    rec_masked(:,:,:,slicenum) = imoverlay(rec(:,:,:,slicenum),~max_mask_topdown,[0 0 0]);
end

%% View a 3d surface
% figure;
% meanrec = mean(rec,4);
% D = imresize3(meanrec,0.25);
scalefactor3d = 0.25;
small_ch1 = (imresize3(squeeze(rec_masked(:,:,1,:)),scalefactor3d));
small_ch2 = (imresize3(squeeze(rec_masked(:,:,2,:)),scalefactor3d));
small_ch3 = (imresize3(squeeze(rec_masked(:,:,3,:)),scalefactor3d));

%%
% % Scatter3
% small_linear = cat(2,small_ch1(:), small_ch2(:), small_ch3(:));
small_linear = small_ch1(:);
msk = find(small_linear(:,1)>50);
pixels = small_linear(msk,:);
pixels = histeq(pixels);
[xx,yy,zz] = ind2sub(size(small_ch1),msk);
scatter3(xx,yy,zz,6,'Marker','*','MarkerEdgeAlpha',0.1,'CData',pixels,'MarkerEdgeColor','flat','MarkerFaceColor','none');
axis vis3d
daspect([1,1,1]);
view(10,10);
% D = imresize3(meanrec,0.25);
% Ds = smooth3(small_ch1);
% clear meanrec
%%
% render = VolumeRender()
% emission = Volume(squeeze(dat(:,:,:,2)));
% render.VolumeEmission=emission;
% render.VolumeAbsorption=emission;
% render.ImageResolution=[size(emission.Data,2), size(emission.Data,1)];
% elementSizeUm = [1;1;1];
% % reflection=Volume([1,1,1;1,1,1;1,1,1]);
% % render.ElementSizeUm(elementSizeUm);
% render.VolumeReflection = emission;
% render.VolumeIllumination = emission;
% render.DistanceToObject=10;
% render.FocalLength=3;
% render.Color=[1,0,0];
% render.VolumeIllumination=Volume(HG(64));
%%
% rendered_image = zeros([size(emission.Data,2), size(emission.Data,1), 3, 360]); 
% nSteps=50;
% angle=360/nSteps;
% for i=1:nSteps
%     render.rotate(0,0,angle);
%     rendered_image(:,:,:,i) = render.render();
% end

%%
body_threshold = 7;
signal_threshold = 30;
clf
smooth_signal = smooth3(small_ch1);
hiso = patch(isosurface(smooth_signal,signal_threshold,smooth_signal),'FaceColor','red','FaceAlpha',1,'EdgeColor','none'); 
isonormals(smooth_signal,hiso);
smoothed_body = smooth3(mean(cat(4,small_ch1, small_ch2),4));
bodypatch = patch(isosurface(smoothed_body,body_threshold,smoothed_body),'FaceAlpha',0.21,'EdgeColor','none'); 
isonormals(smoothed_body, bodypatch);
bodypatch.FaceColor = [0.3 0.3 0.3];

rotate([hiso bodypatch],[1 0 0], 90);
rotate([hiso bodypatch],[0 1 0], 90);
grid on
% c = isocolors(small_ch1,small_ch2,small_ch3, bodypatch);
% h = rgb2hsv(c/255);
% h(:,3) = h(:,3)*16;
% h(:,2) = h(:,2)/2;
% c = hsv2rgb(h);
% bodypatch.FaceVertexCData = c(:,3);

view(30,30);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);

% axis([1,256,1,256,1,256]);
daspect([1,1,1]);
axis vis3d
fixed_axis = axis;
axis(fixed_axis);
lgt = lightangle(30,30);
% lgt = light('Position',[-1 0 0],'Style','infinite');
% lgt = lightangle('Position',[250 0 250])
% camlight
clight = camlight(-50,60,'infinite');
% clight2 = camlight(270,0,'infinite');
lighting gouraud
%% Auto rotate and record video
% rotateby = 2;
% h=waitbar(0,'Rotating');
% for a=1:360/rotateby
%     rotate([hiso bodypatch],[0 1 0], -rotateby);
%     M(a) = getframe(gcf);
%     waitbar(a/360);
% end
% close(h);
% myVideo = VideoWriter([VideoPath,VideoName(1:end-4),'_3dmodel.mp4'],'MPEG-4');
% open(myVideo);
% writeVideo(myVideo,M);
% close(myVideo);