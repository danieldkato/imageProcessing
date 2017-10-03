% demo_ddk.m:

%% DOCUMENTATION TABLE OF CONTENTS:      
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-10-03


%% I. OVERVIEW:
% This script is a wrapper for the NoRMCorre motion-correction algorithm.
% In addition to running the motion correction on an input file, it
% converts and saves the motion-corrected movies as tiffs, saves a .mat
% file including the motion metrics used to evaluate the quality of the
% motion correction, and writes a metadata JSON file including motion
% correction parameters and version information for input files, output
% files, and software dependencies.


%% II. REQUIREMENTS:
% 1) NoRMCorre, available at https://github.com/simonsfoundation/NoRMCorre
% 2) JSONlab, available at https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files
% 3) getLastCommit.m, available at https://github.com/danieldkato/utilities


%% III. INPUTS:
% This script is not a function and has no formal inputs, but requires the
% user to specify a path to a parameters .json file in the code below (see
% comments). For an example of how this file should be formatted, see
% https://github.com/danieldkato/analysis_code/blob/master/motion_correction/normcorre/egparms.json


%% IV. OUTPUT:
% This script is not a function and has no formal return, but saves the
% following to secondary storage:

% 1) a TIFF of the rigid motion-corrected data
% 2) a TIFF of the non-rigid motion-corrected data
% 3) a .mat file containing a struct including NoRMCorre's built-in motion
%    metrics for the raw, rigid motion-corrected, and non-rigid
%    motion-corrected data
% 4) a .JSON metadata file including motion correction parameters and
%    checksums for input files, output files, and software dependencies.


%% TODO:
% This code will have to be modified for longer movies that can't be loaded
% into memory; currently, the code for plotting the motion metrics requires
% that the raw data be loaded into the workspace. 


%% Load parameters:
clear
gcp;

S = loadjson('C:\Users\Dank\Desktop\ncparms.json'); % specify parameters file here

tiffIdx = cell2mat(cellfun(@(x) strcmp(x.input_name,'tiff_to_correct'), S.inputs, 'UniformOutput', false)); 
name = S.inputs{tiffIdx}.path; 
[status, cmdout] = system(['fciv.exe -sha1 ' name]);
S.inputs{tiffIdx}.sha1 = cmdout(end-length(name)-41:end-length(name)-2);

tic; Y = read_file(name); toc; % DDK 2017-09-30: this is fine for small test movies, but will have to be replaced with the name of a .raw file for larger movies
Y = single(Y);                 % convert to single precision 
T = size(Y,ndims(Y));
Y = Y - min(Y(:));


%% set parameters (first try out rigid motion correction)
params = S.params;
options_rigid = NoRMCorreSetParms('d1',params.d1,'d2',params.d2,'bin_width',params.bin_width,'max_shift',params.max_shift,'us_fac',params.us_fac);


%% perform rigid motion correction
tic; [M1,shifts1,template1] = normcorre(Y,options_rigid); toc


%% now try non-rigid motion correction (also in parallel)
options_nonrigid = NoRMCorreSetParms('d1',params.d1,'d2',params.d2,'grid_size',params.grid_size,'mot_uf',params.mot_uf,'bin_width',params.bin_width,'max_shift',params.max_shift,'max_dev',params.max_dev,'us_fac',params.us_fac);
tic; [M2,shifts2,template2] = normcorre_batch(Y,options_nonrigid); toc


%% compute metrics

nnY = quantile(Y(:),0.005); % DDK 2017-09-30: this won't work for larger movies that we can't load into memory; maybe use the mean image instead?
mmY = quantile(Y(:),0.995); % DDK 2017-09-30: this won't work for larger movies that we can't load into memory; maybe use the mean image instead?

[cY,mY,vY] = motion_metrics(Y,10);
[cM1,mM1,vM1] = motion_metrics(M1,10);
[cM2,mM2,vM2] = motion_metrics(M2,10);
T = length(cY);

%% plot metrics
figure;
    ax1 = subplot(2,3,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,3,2); imagesc(mM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    ax3 = subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2,ax3],'xy')
    
    
%% plot shifts        

shifts_r = squeeze(cat(3,shifts1(:).shifts));
shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T); % DDK 2017-09-30: ndims(Y) and T won't be defined for longer movies not loaded into memory
shifts_x = squeeze(shifts_nr(:,1,:))';
shifts_y = squeeze(shifts_nr(:,2,:))';

patch_id = 1:size(shifts_x,2);
str = strtrim(cellstr(int2str(patch_id.')));
str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold') % DDK 2017-09-30: T won't b defined for longer movies that can't be loaded into memory
            set(gca,'Xtick',[]) 
    ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')

    
%% Save output:
t = now;
dstr = datestr(t, 'yyyy-mm-dd');
tstr = datestr(t, 'HH:MM:SS');
dtstr = strrep([dstr '_' tstr],':','-');
dirName = ['mc_output_' dtstr];
mkdir(dirName);
old = cd(dirName);

% Save motion-corrected movies as tiffs:
rmcName = [cd filesep 'rigidMC_' dtstr '.tif'];
saveastiff(M1, rmcName);
nrmcName = [cd filesep 'nonrigidMC_' dtstr '.tif'];
saveastiff(M2, nrmcName);

% Assemble motion metrics into a struct and save as .mat
MM.Uncorrected.CorrCoeffs = cY;
MM.Uncorrected.MeanImg = mY;
MM.Uncorrected.Gradient = vY;
MM.Rigid.CorrCoeffs = cM1;
MM.Rigid.MeanImg = mM1;
MM.Rigid.Gradient = vM1;
MM.Nonrigid.CorrCoeffs = cM2;
MM.Nonrigid.MeanImg = mM2;
MM.Nonrigid.Gradient = vM2;
metricsName = [cd filesep 'motion_metrics_' dtstr '.mat'];
save(metricsName, 'MM');

% Compute and save checksums for output files:
P(1).fieldName = 'rigidMCtiff';
P(1).file = rmcName;
P(2).fieldName = 'nonrigidMCtiff';
P(2).file = nrmcName;
P(3).fieldName = 'motionMetrics';
P(3).file = metricsName;

for i = 1:length(P)
    fname = P(i).file;
    [status, cmdout] = system(['fciv.exe -sha1 ' fname]);
    sha1 = cmdout(end-length(fname)-41:end-length(fname)-2);
    
    S.outputs(i).output_name = P(i).fieldName;
    S.outputs(i).path = fname;
    S.outputs(i).sha1 = sha1;
end

% Add a few miscellaneous fields to metadata structure:
S.date = dstr;
S.time = tstr;
[err, hname] = system('hostname');
S.host_name = strtrim(hname);
S.notes = '';

% Find software dependencies and get version information where available:
ST = dbstack('-completenames');
S.code.main_script.path = ST(1).file;
[warn, latestCommit] = getLastCommit(S.code.main_script.path);
S.code.main_script.lastCommit = latestCommit;
if ~isempty(warn)
    S.code.main_script.warnings = strtrim(warn);
end
Deps = inmem('-completenames');
Deps = Deps(cellfun(@(x) ~contains(x, matlabroot),Deps));
for d = 2:length(Deps)
    clear warn
    S.code.dependencies(d).path = Deps{d};
    [warn, latestCommit] = getLastCommit(Deps{d});
    S.code.dependencies(d).lastCommit = latestCommit;
    if ~isempty(warn)
        S.code.dependencies(d).warnings = strtrim(warn);
    end
end

% Save metadata:
savejson('', S, ['MC_metadata_' dtstr '.json']);
cd(old);
%% plot a movie with the results

%{
figure;
for t = 1:1:T
    subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(122);imagesc(M2(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.02);
end
%}