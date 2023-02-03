%% Parameter setting %%

% % % % % % Confirm your parameters % % % % % % %

ROI = 65; d_ROI = 129; 
ROI_z = 17; d_ROI_z = 30;
% ROI_z = 5; d_ROI_z = 10; 

param.nx = ROI; % number of voxels
param.ny = ROI;
param.nz = ROI_z;

param.sx = ROI*(0.388*460/880); % mm (real size)
param.sy = ROI*(0.388*460/880); % mm
param.sz = ROI_z*(0.388*460/880); % mm

%The real detector panel pixel density (number of pixels)
param.nu = d_ROI;		% number of pixels
param.nv = d_ROI_z;

% Detector setting (real size)
param.su = d_ROI*0.388;	% mm (real size)
param.sv = d_ROI_z*0.388;     % mm

% X-ray source and detector setting
param.DSD = 880;    %  Distance source to detector 
param.DSO = 460;	%  X-ray source to object axis distance

% angle setting
param.dir = -1;   % gantry rotating direction (clock wise/ counter clockwise)
param.dang = 360/600; % angular step size (deg)
param.deg = param.dang:param.dang:360; % you can change
param.deg = param.deg*param.dir;
param.nProj = length(param.deg);

param.parker = 0; % data with 360 deg -> param.parker = 0 , data less than 360 deg -> param.parker=1 

% % % % % % Confirm your parameters % % % % % % %
 
% filter='ram-lak','cosine', 'hamming', 'hann' 
% param.filter='hann'; % high pass filter

param.dx = param.sx/param.nx; % single voxel size
param.dy = param.sy/param.ny;
param.dz = param.sz/param.nz;
param.du = param.su/param.nu;
param.dv = param.sv/param.nv;

param.off_u = param.du*0.25; param.off_v = 0; % detector rotation shift (real size)

% % % Geometry calculation % % %
param.xs = [-(param.nx-1)/2:1:(param.nx-1)/2]*param.dx;
param.ys = [-(param.ny-1)/2:1:(param.ny-1)/2]*param.dy;
param.zs = [-(param.nz-1)/2:1:(param.nz-1)/2]*param.dz;

param.us = (-(param.nu-1)/2:1:(param.nu-1)/2)*param.du + param.off_u;
param.vs = (-(param.nv-1)/2:1:(param.nv-1)/2)*param.dv + param.off_v;

param.interptype = 'linear'; % 'linear', 'nearest'

% % % % % % Confirm your parameters % % % % % % %
% Only for Matlab version above 2013b with parallel computing toolbox: Speed dependent on your GPU
% You don't need to install CUDA, but install latest graphics driver.
% only Nvidia GPU cards can use this. otherwise please "param.gpu=0"
% This option is semi-GPU code using built-in Matlab GPU functions: several times faster than CPU
param.gpu = 0;


















