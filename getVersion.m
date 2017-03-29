function [path, commit] = getVersion(path)
    [pathstr, name, ext] = fileparts(path);
    old = cd(pathstr);
    [status, commit] = system(strcat(['git rev-list --abbrev-commit -1 HEAD ', name, ext]));
    cd(old);
    
    if isempty(commit) || strcmp(commit(1:5),'fatal') == 1
        commit = '(no version info available)';
    else
        commit = strcat(['(commit ', commit, ')']);
    end
end