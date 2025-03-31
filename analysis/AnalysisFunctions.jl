

function readdir_filt(folder_path)


    directories = filter!(e->e!=".DS_Store",readdir(folder_path))
    return directories

end