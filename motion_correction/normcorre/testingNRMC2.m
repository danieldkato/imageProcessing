try

%% get the files for this session

% link to DDK test movie:
% \\10.112.43.46\Public\dank\MultiSens\raw\5036-2\2P\2017-09-25\site1\grab01\2P\5036-2_2017-09-25_site1_grab01_500-2000.tif

%make a big loop
fpath = '/home/georgia/ToBeMotionCorrected/PM04-2/170404/';
cd(fpath)
mouse = 'PM04-2';
filelist = dir([fpath '*.raw']);

%%
for f=1:length(filelist)

%% read big tiff

fname=filelist(f).name;
sesh = strtok(fname,'.');
gcp

% Y = ScanImageTiffReader(fname).data();
% Y = double(Y);

% %% compute metrics for original
% nnY = quantile(Y(:),0.005);
% mmY = quantile(Y(:),0.995);
% 
% %consider adding portion for increasing template init batch if the first
% %section has a weird correlation
% 
% [cY,mY,vY] = motion_metrics(Y,10,10);
% d = ndims(Y);
% T = length(cY);

% %% set parameters (first try out rigid motion correction)
% 
% options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',150,'max_shift',25,'us_fac',50);
% 
% %% perform motion correction
% tic; [M1,shifts1,template1] = normcorre(Y,options_rigid); toc
% 
% %% calculate motion metrics
% [cM1,mM1,vM1] = motion_metrics(M1,10,10);

%% perform non-rigid motion correction

%left off with these for troubleshooting syp/tdt mice
% options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[128,128],...
%     'mot_uf',10,'bin_width',50,'max_shift',25,'max_dev',20,'us_fac',50,...
%     'phase_flag',true,'iter',3,'init_batch',100,...
%     'output_type','memmap','use_parallel',true);

%eftychios params for first test axon video, iter 2 because small file
%(1000 frames)
%options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[128,128],'mot_uf',8,...
 %   'bin_width',50,'max_shift',25,'max_dev',10,'us_fac',50,'phase_flag',true,'iter',2);

% options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[128,128],'mot_uf',8,...
%    'bin_width',50,'max_shift',25,'max_dev',15,'us_fac',50,'phase_flag',true,'iter',1,'output_type','memmap','use_parallel',true);
options_nonrigid = NoRMCorreSetParms('d1',512,'d2',512,'grid_size',[128,128],'mot_uf',8,...
   'bin_width',50,'max_shift',25,'max_dev',15,'us_fac',50,'phase_flag',true,'iter',1,'output_type','memmap','use_parallel',true);

tic; 
% [NR,shifts2,template2] = normcorre_batch(fname,options_nonrigid);
normcorre_batch(fname,options_nonrigid);
toc

    


% %% compute metrics for NR
% [cNR,mNR,vNR] = motion_metrics(NR,10,10);
% 
% %% plot metrics
% figure;
%     ax1 = subplot(2,2,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off;
%         title('mean raw data','fontsize',14,'fontweight','bold')
%     ax2 = subplot(2,2,2); imagesc(mNR,[nnY,mmY]);  axis equal; axis tight; axis off;
%         title('mean raw data','fontsize',14,'fontweight','bold')
%     linkaxes([ax1,ax2],'xy')
%     
%     subplot(2,2,3); plot(1:T,cY,1:T,cNR); legend('raw data','non-rigid'); 
%         title('correlation coefficients','fontsize',14,'fontweight','bold')
%     subplot(2,2,4); scatter(cY,cNR); hold on; plot([0.9*min(cY),1.05*max(cNR)],[0.9*min(cY),1.05*max(cNR)],'--r'); axis square;
%         xlabel('raw data','fontsize',14,'fontweight','bold'); 
%         ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
%     
% %% plot shifts        
% shifts_r = horzcat(shifts1(:).shifts)';
% shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
% shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
% shifts_x = squeeze(shifts_nr(:,1,:))';
% shifts_y = squeeze(shifts_nr(:,2,:))';
% 
% patch_id = 1:size(shifts_x,2);
% str = strtrim(cellstr(int2str(patch_id.')));
% str = cellfun(@(x) ['patch # ',x],str,'un',0);
% 
% figure;
%     ax1 = subplot(311); plot(1:T,cY,1:T,cNR); legend('raw data','non-rigid'); 
%         title('correlation coefficients','fontsize',14,'fontweight','bold')
%             set(gca,'Xtick',[])
%     ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
%             set(gca,'Xtick',[])
%     ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
%             xlabel('timestep','fontsize',14,'fontweight','bold')
%     linkaxes([ax1,ax2,ax3],'x')
% 
%% save a simple tiff for checking in imagej


NR = matfile('motion_corrected.mat');

%select random set of frames to check by eye (just going to say all vid
%have more than 45000 frames for now
randFrame = randi(44000);
saveastiff(uint16(single(NR.Y(:,:,randFrame:(randFrame+999)))),[fpath sesh '_MC-' num2str(randFrame) '.tif']); %,tifoptions)

%% catch any errors and text

%{
load('phone #')
textme(phone,['Motion correction for ' filelist(f).name 'complete!'])

end %end file list loop


catch ME
    % load('phone #')
    textme(phone,[':( Error ' ME.identifier ])
    rethrow(ME)
end
%}