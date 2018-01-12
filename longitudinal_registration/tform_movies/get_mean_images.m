function Movies = get_mean_images(Movies)
%% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-11-29

%% I. OVERVIEW:
% This function goes through each element of the struct array Movies, each
% of which represents a movie to be concentenated by
% register_longitudinal_data, and computes its mean image. This will be
% useful for computing what transformations to use to register all of the
% movies to a reference movie (which will be the mean image of the first
% movie).


%% II. REQUIREMENTS:
% 1) The MATLAB Parallel Computing toolbox


%% III. INPUTS:
% 1) Movies - 1 x m struct array where each element represents a movie to be
%    concatenated. For more detail on how this should be formatted, see the
%    OUTPUTS section of the documentation for initialize_movie_struct.m.


%% IV. OUTPUTS:
% 1) Movies - same as the input, but now each element also includes a
%    non-empty 'mean_img' field, an m-by-n matrix of image data specifying
%    the mean image of the corresponding movie, where m is the movie width
%    in pixels and n is the movie height in pixels.


%%
disp('Finding average image for each movie...');
parfor m = 1:length(Movies)
    disp(['finding average image for movie ' num2str(m)]);
    var_name = Movies(m).var_name;
    Movies(m).mean_img = mean(Movies(m).matfile.(var_name),3);
    disp(['done finding average image for movie ' num2str(m)]);
end
disp('... done finding average image for each movie.');

end