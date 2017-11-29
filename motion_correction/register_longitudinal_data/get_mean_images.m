function Movies = get_mean_images(Movies)

disp('Finding average image for each movie...');
parfor m = 1:length(Movies)
    disp(['finding average image for movie ' num2str(m)]);
    var_name = Movies(m).var_name;
    Movies(m).mean_img = mean(Movies(m).matfile.(var_name),3);
    disp(['done finding average image for movie ' num2str(m)]);
end
disp('... done finding average image for each movie.');

end