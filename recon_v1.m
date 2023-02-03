clear;
close all;
clc

warning('off','all')

NoiseOn = 0;    % 0- without noise, 1- with noise



ROIx = 129; ROIy = 129; tt = 64; ROIz = 20; %ROIz = 33; 
pixel_size = 0.2028;
lesioncoef = 0.238;
Ninput = 6914;
diameter = 1;
enoise_var = 0; % if there is electronic noise
elipseA = 1;
elipseB = 1;

VGF = 50;
startm = 0;
iterateM = 2;
Nimg =5000;
savetag = 'Central_caseM';
dataRate = 1;


%% 3D image
W = zeros(ROIx,ROIy,ROIz);
for xx = 1:ROIx
    for yy = 1:ROIy
        for zz = 1:ROIz
            if sqrt((xx-floor((ROIx)/2+1))^2+(yy-floor((ROIy)/2+1))^2+(zz-floor((ROIz)/2+1))^2)<=tt
                W(xx,yy,zz) = 1;
            end
        end
    end
end

ROI_filter = 257;

delta_f_filter = 1/(ROI_filter*pixel_size);
kf = -floor(ROI_filter/2) : floor(ROI_filter/2);
lf = -floor(ROI_filter/2) : floor(ROI_filter/2);
mf = -floor(ROI_filter/2) : floor(ROI_filter/2);


f_filter = zeros(ROI_filter,ROI_filter,ROI_filter);
for x=1:length(kf)
    for y=1:length(lf)
        for z=1:length(mf)
            f_filter(x,y,z) = delta_f_filter*sqrt(kf(x)^2+lf(y)^2+mf(z)^2);
        end
    end
end

beta = 3;
filter = 1./(f_filter.^(beta/2));
filter(floor(ROI_filter/2)+1,floor(ROI_filter/2)+1,floor(ROI_filter/2)+1) = 0;
filter(floor(ROI_filter/2)+1,floor(ROI_filter/2)+1,floor(ROI_filter/2)+1) = 2*max(max(max(filter)));


%% Image generation
for m = 1:iterateM
    hann_G1  = zeros(floor(ROIx/2)+1,floor(ROIy/2)+1,17,Nimg);
    ramp_G1  = zeros(floor(ROIx/2)+1,floor(ROIy/2)+1,17,Nimg);
    
    
    time = 0;
    for ii=1:Nimg
        tic
        %         matlabpool open;
        bg_temp = ifftn(ifftshift(fftshift(fftn(randn(ROI_filter,ROI_filter,ROI_filter))).*filter));
        img_temp = bg_temp(floor(ROI_filter/2)+1-tt:floor(ROI_filter/2)+1+tt,floor(ROI_filter/2)+1-tt:floor(ROI_filter/2)+1+tt,floor(ROI_filter/2)+1-tt:floor(ROI_filter/2)+1+tt);
        temp = img_temp;
        sorted = sort(temp(:),'descend');
        threshold_val = sorted(round(length(sorted)*(VGF/100)));
        img = zeros(ROIx,ROIy,ROIz);
        for x=1:ROIx
            for y=1:ROIy
                for z=1:ROIz
                    if temp(x,y,z) > threshold_val
                        img(x,y,z) = 0.233;
                    else
                        img(x,y,z) = 0.194;
                    end
                end
            end
        end
        
        for x=1:ROIx
            for y=1:ROIy
                for z=1:ROIz
                    if sqrt((x-floor(ROIx/2+1))^2+(y-floor(ROIy/2+1))^2+(z-floor(ROIz/2+1))^2) < diameter/pixel_size/2
                        img(x,y,z) = lesioncoef;
                    end
                end
            end
        end
        
        img = img.*W;
        
        % 2. Projection
        
        ParamSetting;
        proj = zeros(param.nu, param.nv, param.nProj,'single');
        parfor k = 1:param.nProj
            proj(:,:,k) = projection(img,param,k);
        end
        
%         ParamSetting;
        proj_quantum = zeros(param.nu,param.nv,param.nProj);
%         N0=6914;
        N0=Ninput;
        Nin=N0*exp(zeros(size(proj_quantum)));
        proj_quantum=10*(log(Nin)-log(poissrnd(N0*ones(size(proj_quantum))))); % x10 because our unit is mm but unit of mu is cm
        
        
        proj = proj + proj_quantum+ sqrt(enoise_var) * randn(size(proj));
        
        % 3. Filtered projection
        
        param.filter='hann';
        proj_filtered_hann  = filtering(proj,param); % vertical: y, horizontal: z
        param.filter='ram-lak';
        proj_filtered_ramp  = filtering(proj,param);
        
        % 4. Sinc Interpolation
        
        % 5. Reconstruction
        ParamSetting65recon;
        Recon_hann_linear = 0; Recon_shepp_linear = 0; Recon_ramp_linear = 0;
        [xx,yy] = meshgrid(param.xs,param.ys);
        
        parfor iview = 1:param.nProj
            angle_rad = param.deg(iview)/360*2*pi;
            
            
            rx = xx.*cos(angle_rad-pi/2) + yy.*sin(angle_rad-pi/2);
            ry = -xx.*sin(angle_rad-pi/2) + yy.*cos(angle_rad-pi/2);
            
            pu = single(((rx.*(param.DSD)./(ry + param.DSO))+param.us(1))/(-param.du) + 1);
            Ratio = (single(param.DSO.^2./(param.DSO+ry).^2));
            
            vol_hann_linear  = zeros(param.nx,param.ny,param.nz,'single');
            vol_ramp_linear  = zeros(param.nx,param.ny,param.nz,'single');
            
            
            for iz = 1:param.nz
                pv = single(((param.zs(iz)*(param.DSD)./(ry + param.DSO))-param.vs(1))/param.dv+1);
                vol_hann_linear(:,:,iz)  = (Ratio.*interp2(proj_filtered_hann(:,:,iview)',pu,pv,param.interptype));
                vol_ramp_linear(:,:,iz)  = (Ratio.*interp2(proj_filtered_ramp(:,:,iview)',pu,pv,param.interptype));
            end
            
            vol_hann_linear(isnan(vol_hann_linear))=0;
            vol_ramp_linear(isnan(vol_ramp_linear))=0;
            
            
            Recon_hann_linear = Recon_hann_linear + vol_hann_linear;
            Recon_ramp_linear = Recon_ramp_linear + vol_ramp_linear;
            
            
            
            
        end
        
        hann_G1(:,:,:,ii) = Recon_hann_linear(:,:,:);
        ramp_G1(:,:,:,ii) = Recon_ramp_linear(:,:,:);
        
        
        clc
        %         matlabpool close;
        temptime = toc
        time = time + temptime;
        remaintime = time/ii * (Nimg*(1+dataRate)-ii);
        fprintf('G1 생성중 ii = %d, 예상 남은시간은 %2d일 %3d시간 %3d분 %3.2f초입니다.',...
            ii, floor(remaintime/3600/24),floor(mod(remaintime/3600,24)),floor(mod(remaintime,3600)/60),mod(remaintime,60));
        
    end
    
    % g0 set generation
    
    hann_G0  = zeros(floor(ROIx/2)+1,floor(ROIy/2)+1,floor(ROIz/2)+1,Nimg*dataRate);
    ramp_G0  = zeros(floor(ROIx/2)+1,floor(ROIy/2)+1,floor(ROIz/2)+1,Nimg*dataRate);
    
    
    time=0;
    for ii=1:Nimg*dataRate
        tic
        %         matlabpool open;
        bg_temp = ifftn(ifftshift(fftshift(fftn(randn(ROI_filter,ROI_filter,ROI_filter))).*filter));
        img_temp = bg_temp(floor(ROI_filter/2)+1-tt:floor(ROI_filter/2)+1+tt,floor(ROI_filter/2)+1-tt:floor(ROI_filter/2)+1+tt,floor(ROI_filter/2)+1-tt:floor(ROI_filter/2)+1+tt);
        temp = img_temp;
        sorted = sort(temp(:),'descend');
        threshold_val = sorted(round(length(sorted)*(VGF/100)));
        img = zeros(ROIx,ROIy,ROIz);
        for x=1:ROIx
            for y=1:ROIy
                for z=1:ROIz
                    if temp(x,y,z) > threshold_val
                        img(x,y,z) = 0.233;
                    else
                        img(x,y,z) = 0.194;
                    end
                end
            end
        end
        
        img = img.*W;
        
        % 2. Projection
        
        ParamSetting;
        proj = zeros(param.nu, param.nv, param.nProj,'single');
        parfor i = 1:param.nProj
            proj(:,:,i) = projection(img,param,i);
        end
        
        ParamSetting;
        proj_quantum = zeros(param.nu,param.nv,param.nProj);
        N0=Ninput;
        Nin=N0*exp(zeros(size(proj_quantum)));
        proj_quantum=10*(log(Nin)-log(poissrnd(N0*ones(size(proj_quantum))))); % x10 because our unit is mm but unit of mu is cm
        
        
        proj = proj + proj_quantum+ sqrt(enoise_var) * randn(size(proj));
        
        % 3. Filtered projection
        
        param.filter='hann';
        proj_filtered_hann  = filtering(proj,param); % vertical: y, horizontal: z
        param.filter='ram-lak';
        proj_filtered_ramp  = filtering(proj,param);
        
        % 4. Sinc Interpolation
        
        % 5. Reconstruction
        ParamSetting65recon;
        Recon_hann_linear = 0; Recon_shepp_linear = 0; Recon_ramp_linear = 0;
        [xx,yy] = meshgrid(param.xs,param.ys);
        
        parfor iview = 1:param.nProj
            angle_rad = param.deg(iview)/360*2*pi;
            
            
            rx = xx.*cos(angle_rad-pi/2) + yy.*sin(angle_rad-pi/2);
            ry = -xx.*sin(angle_rad-pi/2) + yy.*cos(angle_rad-pi/2);
            
            pu = single(((rx.*(param.DSD)./(ry + param.DSO))+param.us(1))/(-param.du) + 1);
            Ratio = (single(param.DSO.^2./(param.DSO+ry).^2));
            
            vol_hann_linear  = zeros(param.nx,param.ny,param.nz,'single');
            vol_ramp_linear  = zeros(param.nx,param.ny,param.nz,'single');
            
            
            for iz = 1:param.nz
                pv = single(((param.zs(iz)*(param.DSD)./(ry + param.DSO))-param.vs(1))/param.dv+1);
                vol_hann_linear(:,:,iz)  = (Ratio.*interp2(proj_filtered_hann(:,:,iview)',pu,pv,param.interptype));
                vol_ramp_linear(:,:,iz)  = (Ratio.*interp2(proj_filtered_ramp(:,:,iview)',pu,pv,param.interptype));
            end
            
            vol_hann_linear(isnan(vol_hann_linear))=0;
            vol_ramp_linear(isnan(vol_ramp_linear))=0;
            
            
            Recon_hann_linear = Recon_hann_linear + vol_hann_linear;
            Recon_ramp_linear = Recon_ramp_linear + vol_ramp_linear;
            
            
            
            
        end
        
        hann_G0(:,:,:,ii) = Recon_hann_linear(:,:,:);
        ramp_G0(:,:,:,ii) = Recon_ramp_linear(:,:,:);
        
        
        clc
        temptime = toc
        time = time + temptime;
        remaintime = time/ii * (Nimg*dataRate-ii);
        fprintf('G0 생성중 ii = %d, 예상 남은시간은 %2d일 %3d시간 %3d분 %3.2f초입니다.\n',...
            ii, floor(remaintime/3600/24),floor(mod(remaintime/3600,24)),floor(mod(remaintime,3600)/60),mod(remaintime,60));
        
    end
    save(['G:\CT\ObserverData\Observer_3D',savetag, int2str(startm+m), '.mat'], 'hann_G0', 'hann_G1', 'ramp_G0', 'ramp_G1','-v7.3')
end

