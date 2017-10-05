% ddk_normcorre_wrapper.m:

%% DOCUMENTATION TABLE OF CONTENTS:      
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-10-05


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
% https://github.com/danieldkato/analysis_code/blob/master/motion_correction/normcorre/mc_params.json


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

% add NoRMCorre and utilities to path; hopefully this won't be necessary
% after I figure out issues with pathdef.m
p = genpath('~/Documents/MATLAB/NoRMCorre');
addpath(p);
p2 = genpath('/mnt/nas2/homes/dan/libraries/utilities');
addpath(p2);
p3 = genpath('~/Documents/MATLAB/Add-ons/Toolboxes/manur-MATLAB-git-72c72b9');
addpath(p3);

S = loadjson('/mnt/nas2/homes/dan/libraries/image_processing/motion_correction/normcorre/mc_params.json'); % specify parameters file here

tiffIdx = cell2mat(cellfun(@(x) strcmp(x.input_name,'tiff_to_correct'), S.inputs, 'UniformOutput', false)); 
name = S.inputs{tiffIdx}.path; 
if ispc
    [status, cmdout] = system(['fciv.exe -sha1 ' name]);
elseif isunix
    [status, cmdout] = system(['sha1sum ' name]);
end
S.inputs{tiffIdx}.sha1 = cmdout(end-length(name)-41:end-length(name)-2);

cd(S.outputs.output_directory)
S = rmfield(S,'outputs');

tic; Y = read_file(name); toc; % DDK 2017-09-30: this is fine for small test movies, but will have to be replaced with the name of a .raw file for larger movies
Y = single(Y);                 % convert to single precision 
T = size(Y,ndims(Y));
Y = Y - min(Y(:));


%% set parameters (first try out rigid motion correction)

% Will read parameter names and values into the string that you would type
% into the MATLAB command line to invoke NoRMCorre; this string will be
% passed to eval() to create the options structure. (The reason for this is
% that it allows the user to add or remove parameters from the parameters
% JSON file without having to change the code).

nr_param_names = fieldnames(S.params);
r_param_names = {'d1','d2','bin_width','max_shift','us_fac'};

set_r_params_str = 'NoRMCorreSetParms(';
set_nr_params_str = set_r_params_str;

% go through every parameters in the parameter file:
rigid_param_ctr = 0;
for s = 1:length(nr_param_names)
    
    name = nr_param_names{s}; % get its name
    val = S.params.(nr_param_names{s}); % get its value
    
    % convert the value into a string:
    if isnumeric(val) && length(val) == 1
        vStr = num2str(val);
    elseif isnumeric(val) && length(val) > 1
        vStr = '[';
        for v = 1:length(val)
            vStr = [vStr num2str(val(v)) ' '];
        end
        vStr = [vStr ']'];
    end    
    
    % concat name and value strings:
    kv = ['''' name ''',' vStr];
    if s<length(nr_param_names)
       kv = [kv ',']; 
    end
        
    % add it to the non-rigid mc parameters
    set_nr_params_str = [set_nr_params_str kv];
    
    % check if it is also a rigid mc parameter
    match = cellfun(@(x) strcmp(x, name), r_param_names);
    if sum(match) > 0 
        rigid_param_ctr = rigid_param_ctr+1;
        if rigid_param_ctr == length(r_param_names)
            kv = kv(1:end-1);
        end
        set_r_params_str = [set_r_params_str kv];
    end
    
end
set_r_params_str = [set_r_params_str ')'];
disp(set_r_params_str);
set_nr_params_str = [set_nr_params_str ')'];
disp(set_nr_params_str);

options_rigid = eval(set_r_params_str);


%% perform rigid motion correction
tic; [M1,shifts1,template1] = normcorre(Y,options_rigid); toc


%% now try non-rigid motion correction (also in parallel)
options_nonrigid = eval(set_nr_params_str);
tic; [M2,shifts2,template2] = normcorre_batch(Y,options_nonrigid); toc


%% compute metrics
nnY = quantile(Y(:),0.005); % DDK 2017-09-30: this won't work for larger movies that we can't load into memory; maybe use the mean image instead?
mmY = quantile(Y(:),0.995); % DDK 2017-09-30: this won't work for larger movies that we can't load into memory; maybe use the mean image instead?

[cY,mY,vY] = motion_metrics(Y,10);
[cM1,mM1,vM1] = motion_metrics(M1,10);
[cM2,mM2,vM2] = motion_metrics(M2,10);
T = length(cY);


%% plot metrics
f1 = figure;
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

f2 = figure;
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

% Save figures:
fig1name = [cd filesep 'motion_metrics_fig1' dtstr '.fig'];
fig2name = [cd filesep 'motion_metrics_fig2' dtstr '.fig'];
savefig(f1, fig1name);
savefig(f2, fig2name);

% Compute and save checksums for output files:
P(1).fieldName = 'rigid_mc_tiff';
P(1).file = rmcName;
P(2).fieldName = 'nonrigid_mc_tiff';
P(2).file = nrmcName;
P(3).fieldName = 'motion_metrics';
P(3).file = metricsName;
P(4).fieldName = 'motion_metrics_fig1';
P(4).file = fig1name;
P(5).fieldName = 'motion_metrics_fig2';
P(5).file = fig2name;

disp('Computing output file checksums...');
for i = 1:length(P)
    fname = P(i).file;
    
    tic;
    if ispc
        [status, cmdout] = system(['fciv.exe -sha1 ' fname]);
    elseif isunix
        [status, cmdout] = system(['sha1sum ' fname]);
    end
    
    sha1 = cmdout(end-length(fname)-41:end-length(fname)-2);
    disp([num2str(i) ' out of ' num2str(length(P)) ' checksums complete.']);
    toc;
    
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
disp('Retrieiving software dependencies...')
ST = dbstack('-completenames');
S.code.main_script.path = ST(1).file;
[warn, latestCommit] = getLastCommit(S.code.main_script.path);
S.code.main_script.lastCommit = latestCommit;
if ~isempty(warn)
    S.code.main_script.warnings = strtrim(warn);
end
Deps = inmem('-completenames');
Deps = Deps(cellfun(@(x) ~contains(x,matlabroot), Deps)); % filter out core MATLAB functions - there would be way too many (>600!)
Deps = Deps(cellfun(@(x) ~contains(x,'pathdef.m') & ~contains(x,mfilename), Deps)); % filter out this script and pathdef.m
for d = 1:length(Deps)
    clear warn
    S.code.dependencies(d).path = Deps{d};
    [warn, latestCommit] = getLastCommit(Deps{d});
    S.code.dependencies(d).lastCommit = latestCommit;
    if ~isempty(warn)
        S.code.dependencies(d).warnings = strtrim(warn);
    end
end
disp('... complete.');

disp('Saving metadata...');

% Save metadata:
savejson('', S, ['MC_metadata_' dtstr '.json']);
disp('... complete.');

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