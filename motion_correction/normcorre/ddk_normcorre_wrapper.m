function ddk_normcorre_wrapper(param_file_path)
%% DOCUMENTATION TABLE OF CONTENTS:      
% I. OVERVIEW
% II. USAGE
% III. REQUIREMENTS
% IV. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-18-06


%% I. OVERVIEW:
% This function is a wrapper for the Simons Foundations' NoRMCorre
% motion-correction package. In addition to running motion correction on an
% input movie, it saves to secondary storage the motion-corrected movie in
% a format specified by the user (see NoRMCorre documentation for options),
% a .mat file of the motion metrics used to evaluate the quality of the
% motion correction, and a metadata JSON file including motion correction
% parameters and version information for input files, output files, and
% software dependencies. 


%% II. USAGE:

% In addition to invoking this function in another MATLAB script or from
% the MATLAB command line, it is possible to invoke this function from the
% LINUX command line with the following:

% matlab -nosplash -nodesktop -r "ddk_normcorre_wrapper </path/to/params>"

% where <path/to/params> stands for the path to theinput parameters file.


%% III. REQUIREMENTS:
% 1) The MATLAB package NoRMCorre, available at https://github.com/simonsfoundation/NoRMCorre

%    *** IMPORTANT NOTE ***: if trying to read in a TIFF and write output
%    to a memory-mapped .mat file, it will be necessary to either use the
%    `tiff_2_memmap_mat` git branch of DDK's local clone of the NoRMCorre 
%    repo (available at /mnt/nas2/homes/dan/code_libraries/NoRMCorre) or 
%    apply certain bug fixes to normcorre.m and normcorre_batch.m. As is,
%    NoRMCorre crashes when trying to read in a TIFF and write to a memory-mapped
%    .mat file. For a description of the bug, see commit 187dfb0 on
%    `tiff_2_memmap_mat` branch of DDK's local clone of the NoRMCorre repo.
%    This is ok for now, but I should open an issue on the main NoRMCorre
%    github repository. 

% 2) The MATLAB toolbox JSONlab, available at https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files
% 3) The MATLAB function write_metadata.m, available at https://github.com/danieldkato/utilities
% 4) The MATLAB function generate_mc_dir_name.m, availalbe at https://github.com/danieldkato/imageProcessing/blob/master/motion_correction/generate_mc_dir_name.m


%% IV. INPUTS:
% 1) param_file_path - path to a parameters .json file. Example contents of
%    such a file might be as follows:

% {
% "description":"Motion correction of raw imaging data using NoRMCorre non-rigid motion correction algorithm.",
% "inputs":[
%		{
%		"input_name":"tiff_to_correct",
%		"path":"/mnt/nas2/homes/dan/MultiSens/data/test_movies/5036-2_short/5036-2_short.tif"
%		}
%	],
% "params":{
%	"do_nonrigid":"true",
%    "do_rigid":"false",
%    "save_tiff":"false",
%    "normcorre_params":{
%        "d1":512,
%       "d2":512,
%       "grid_size":[32,32],
%       "mot_uf":4,
%       "bin_width":50,
%       "max_shift":15,
%       "max_dev":3,
%       "us_fac":50,
%       "output_type":"hdf5"
%        }
%	},
%}

% Note that at the moment, this script assumes that the input file is 1)
% large, and therefore doesn't try to load it into memory, and 2) a TIFF.
% Because of this, this script omits certain parts of the original
% NoRMCorre demo, like computing and plotting the motion metrics for the
% raw movie, which requires that the raw movie be small enough to load into
% memory or saved as a .mat. 

% Also, note that if the input file is very large, then the `output_type`
% parameter MUST be set to any value other that `mat` - otherwise,
% NoRMCorre will by default try to save the output as an in-memory variable
% (which the output will be too large for).

% For another example of how this file should be formatted, see
% https://github.com/danieldkato/analysis_code/blob/master/motion_correction/normcorre/mc_params.json


%% V. OUTPUT:
% This function and has no formal return, but saves the following to
% secondary storage:

% 1) A file of the rigid motion-corrected data. The format of this output
%    file is determined by the `output_type` field of the parameters file.

% 2) A .mat file of motion metrics output (see NoRMCorre
%    documentation for more detail on motion metrics).

% Optionally, depending on the user-specified parameters and the format of
% the inputs and outputs, this fucntion may also save:

% 3) A TIFF of the motion-corrected movies. This is specified by the
%    `save_tiff` field of the params file. If set to "true", then even if
%    "output_type" is not set to "tif", then this script will *additionally*
%    save the output as a TIFF. Note, however, that this is extremely
%    time-consuming and not necessary for all subsequent analyses, so it's
%    probably best to keep "false". If "output_type" is set to "tif", then
%    this option has no effect.

% 4) Figures of the motion metrics. This will only be plotted if the input
%    movie is a .mat file. The reason for this is that Eftykios' plots require
%    running motion_metrics() on the input movie, and this seems to fail if
%    the input movie is a TIFF.


%% TODO:
% 1) Would be nice to have some way to also save output as a TIFF for
% quick visual inspection


%% Load parameters:
gcp;

start = tic;
disp('Loading parameters... ');
S = loadjson(param_file_path); % specify parameters file here

% Get name of movie to be motion-corrected (necessary for larger movies
% that can't be loaded into memory):
tiffIdx = cell2mat(cellfun(@(x) strcmp(x.input_name,'tiff_to_correct'), S.inputs, 'UniformOutput', false)); 
input_path = S.inputs{tiffIdx}.path; 
[input_dir input_name input_type] = fileparts(input_path);

% Determine whether a directory called 'motion_correction' exists in the
% directory where the input movie is saved, and if not, create it:
cd(input_dir);
ls = dir();
is_mc = arrayfun(@(x) strcmp(x.name, 'motion_correction'), ls);
if isempty(find(is_mc,1))
    mkdir('motion_correction');
end
orig_dir = cd('motion_correction');

% Create a directory for this specific run of NoRMCorre; need to create
% this BEFORE running normcorre_batch() so memory-mapped output files can
% be saved to correct location:
dirName = generate_mc_dir_name();
disp(cd);
disp(dirName);
mkdir(dirName);
cd(dirName);
base = [cd filesep];

% Initialize metadata struct:
Metadata.inputs(1).path = input_path;
Metadata.params = S.params;
Metadata.outputs = [];

% Determine whether output should be saved as a tiff in addition to any
% other format:
save_tiff = false;
save_tiff_char = S.params.save_tiff;
if strcmp(save_tiff_char, 'true')
    save_tiff = true;
end

disp('... done loading parameters.');

% Code for loading input movie into memory; commenting out because this
% will not work for large movies
%{
% Load image data:
tic; Y = read_file(name); toc; % DDK 2017-09-30: this is fine for small test movies, but will have to be replaced with the name of a .raw file for larger movies
Y = single(Y);                 % convert to single precision 
T = size(Y,ndims(Y));
Y = Y - min(Y(:));
%}

% Get the appropriate output file extension based on the output_type option:
if strcmp(S.params.normcorre_params.output_type, 'tiff')
    ext = 'tif';
    output_name_option = 'tiff_filename';
else
    if strcmp(S.params.normcorre_params.output_type, 'hdf5')
        ext = 'h5';
    else
        ext = S.params.normcorre_params.output_type;
    end
    output_name_option = [ext '_filename'];
end


%% set parameters (first try out rigid motion correction)

% Will combine parameter names and values into a string that one would
% enter into the MATLAB command line to invoke NoRMCorreSetParms; this
% string will be passed to eval() to create the options structure. (The
% reason for this is that it allows the user to add or remove parameters
% from the parameters JSON file without having to change the code).

nr_param_names = fieldnames(S.params.normcorre_params);
disp('nr_param_names:');
disp(nr_param_names);
r_param_names = {'d1','d2','bin_width','max_shift','us_fac','output_type'};

set_r_params_str = 'NoRMCorreSetParms(';
set_nr_params_str = set_r_params_str;

% Go through every parameter in the parameter file:
rigid_param_ctr = 0;
for s = 1:length(nr_param_names)
    
    name = nr_param_names{s}; % get parameter name
    val = S.params.normcorre_params.(name); % get parameter value
    
    % convert the value into a string:
    if isnumeric(val) && length(val) == 1
        vStr = num2str(val);
    elseif isnumeric(val) && length(val) > 1
        vStr = '[';
        for v = 1:length(val)
            vStr = [vStr num2str(val(v)) ' '];
        end
        vStr = [vStr ']'];
    elseif ischar(val)
        vStr = [ '''' val ''''];
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
        if rigid_param_ctr == length(r_param_names) && kv(end) == ','
            kv = kv(1:end-1);
        end
        set_r_params_str = [set_r_params_str kv];
    end
    
end
set_r_params_str = [set_r_params_str ')'];
disp(set_r_params_str);
set_nr_params_str = [set_nr_params_str ')'];
disp(set_nr_params_str);


%% Perform rigid motion correction if specified by user:

% Determine whether user wants to do rigid motion correction:
do_rigid = false;
do_rigid_char = S.params.do_rigid;
if strcmp(do_rigid_char, 'true')
    do_rigid = true;
end

% If requested by user, perform rigid motion correction:
if do_rigid
    
    disp('Performing rigid motion correction...');
    
    rmc_name = [base 'rigidMC.' ext];
    
    % Create the options object:
    options_rigid = eval(set_r_params_str);
    options_rigid.(output_name_option) = rmc_name;
    
    % Do the motion correction using the options object we just created:
    tic; [M1,shifts1,template1] = normcorre(input_path,options_rigid); toc 
    
    % Compute motion metrics:
    if strcmp(S.params.normcorre_params.output_type, 'mat') || strcmp(S.params.normcorre_params.output_type, 'memmap')
        [cM1,mM1,vM1] = motion_metrics(M1,10); 

        % Save the motion metrics to a struct:
        MM.Rigid.CorrCoeffs = cM1;
        MM.Rigid.MeanImg = mM1;
        MM.Rigid.Gradient = vM1;
    end
    
    % Append names of outputs to Metadata struct:
    Metadata.outputs(end+1).path = rmc_name;
    
    % If specified by the user, save the output as a TIFF:
    if save_tiff 
        try
            rmcTifName = [cd filesep 'rigidMC_' dtstr '.tif'];
            saveastiff(M1, rmcTifName);
            Metadata.outputs(end+1).path = rmcTifName;
        catch
            warning('Specified output type incompatible with saveastiff(); will not save output as TIFF.');
        end
    end
    
    disp('... done performing rigid motion correction.');
    
end


%% Now try non-rigid motion correction (also in parallel) if specified by the user:

% Determine whether user has requested non-rigid motion correction:
do_nonrigid = false;
do_nonrigid_char = S.params.do_nonrigid;
if strcmp(do_nonrigid_char, 'true')
    do_nonrigid = true;
end

% If requested by user, perform non-rigid motion correction:
if do_nonrigid
    
    disp('Performing non-rigid motion correction...');
    
    nrmc_name = [base 'nonrigidMC.' ext];
    
    % Create the options object:
    options_nonrigid = eval(set_nr_params_str);
    options_nonrigid.(output_name_option) = nrmc_name;
    
    % Do the motion correction using the options object we just created:
    tic; [M2,shifts2,template2] = normcorre_batch(input_path,options_nonrigid); toc % do the motion correction
    
    % Compute motion metrics:
    if strcmp(S.params.normcorre_params.output_type, 'mat') || strcmp(S.params.normcorre_params.output_type, 'memmap')
        [cM2,mM2,vM2] = motion_metrics(M2,10); 

        %Save motion metrics to a struct:
        MM.Nonrigid.CorrCoeffs = cM2;
        MM.Nonrigid.MeanImg = mM2;
        MM.Nonrigid.Gradient = vM2;
    end
    
    
    % Append names of outputs to Metadata struct:
    Metadata.outputs(end+1).path = nrmc_name;
    
    % If specified by the user, save the output as a TIFF:
    if save_tiff 
        try
            nrmcTifName = [cd filesep 'nonrigidMC_' dtstr '.tif'];
            saveastiff(M2, nrmcTifName);
            Metadata.outputs(end+1).path = nrmcTifName;
        catch
            warning('Specified output type incompatible with saveastiff(); will not save output as TIFF.');
        end
    end
    
    disp(' ... done performing non-rigid motion correction.');
    
end


%% Plot metrics *only if input type is .mat*. Otherwise, running motion_metrics() on the input movie will fail, and this is necessary for Eftykios' plots:

if strcmp(input_type, '.mat')

    disp('Plotting motion metrics...');
    
    [cY,mY,vY] = motion_metrics(Y,10);
    T = length(cY);
    nnY = quantile(M2.Y(:,:,:),0.005); % DDK 2017-09-30: this won't work for larger movies that we can't load into memory; maybe use the mean image instead?
    mmY = quantile(M2.Y(:,:,:),0.995); % DDK 2017-09-30: this won't work for larger movies that we can't load into memory; maybe use the mean image instead?


    f1 = figure;
        ax1 = subplot(2,3,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
        ax2 = subplot(2,2,2); imagesc(mM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
        ax3 = subplot(2,2,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
        subplot(2,2,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
        subplot(2,2,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
            xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
        subplot(2,2,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
            xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
        linkaxes([ax1,ax2,ax3],'xy')


    %plot shifts        
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
    shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
    shifts_nr = reshape(shifts_nr,[],ndims(M2)-1,T); % DDK 2017-09-30: ndims(Y) and T won't be defined for longer movies not loaded into memory
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
        
    % Save figures:
    fig1name = [cd filesep 'motion_metrics_fig1' dtstr '.fig'];
    fig2name = [cd filesep 'motion_metrics_fig2' dtstr '.fig'];
    savefig(f1, fig1name);
    savefig(f2, fig2name);
    
    % Append figures to Metadata:
    Metadata.outputs(end+1).path = fig1name;
    Metadata.outputs(end+1).path = fig2name;

    disp('... done plotting motion metrics.');
    
end

    
%% Save output:

% Save motion metrics as .mat
if strcmp(S.params.normcorre_params.output_type, 'mat') || strcmp(S.params.normcorre_params.output_type, 'memmap')
    metricsName = [cd filesep 'motion_metrics.mat'];
    save(metricsName, 'MM');
    Metadata.outputs(end+1).path = metricsName;
end 

%% Save metadata:

disp('Saving metadata...');
Metadata.duration = [num2str(toc(start)) ' seconds'];
write_metadata(Metadata, 'MC_metadata.json');
disp('... done saving metadata.');

cd(orig_dir);


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