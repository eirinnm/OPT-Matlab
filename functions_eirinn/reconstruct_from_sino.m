function s = reconstruct_from_sino(sino, angles, gpuAvailable)
    % create geometries and projector
    image_width = size(sino,2);
    proj_geom = astra_create_proj_geom('parallel', 1.0, image_width, angles);
    vol_geom = astra_create_vol_geom(image_width,image_width);
    recon_id = astra_mex_data2d('create', '-vol', vol_geom);
    sinogram_id = astra_mex_data2d('create','-sino', proj_geom);

    if gpuAvailable == 1
        cfg = astra_struct('FBP_CUDA');
        proj_id = astra_create_projector('cuda', proj_geom, vol_geom); 
    else
        cfg = astra_struct('FBP');
        proj_id = astra_create_projector('linear', proj_geom, vol_geom); % strip or linear?
    end
%     if exist('sino_mask','var')
%         mask_geom = astra_create_vol_geom(size(sino,1),size(sino,2));
%         disp(mask_geom)
%         mask_id = astra_mex_data2d('create', '-vol', mask_geom);
%         astra_mex_data2d('set',mask_id,sino_mask);
%         cfg.option.SinogramMaskID = mask_id;
%     end
    cfg.ProjectorId = proj_id;
    cfg.ProjectionDataId = sinogram_id;
    cfg.ReconstructionDataId = recon_id;
    
    %cfg.FilterType = 'Ram-Lak';
    %cfg.FilterType = 'shepp-logan';
    astra_mex_data2d('set',sinogram_id,single(sino));
    % run the algorithm
    fbp_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('run', fbp_id);
    s = astra_mex_data2d('get', recon_id);
    % cleanup
    astra_mex_data2d('delete', sinogram_id, recon_id);
    astra_mex_projector('delete', proj_id);
    astra_mex_algorithm('delete', fbp_id);
end


