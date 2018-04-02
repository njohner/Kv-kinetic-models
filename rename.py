import os
import subprocess

basepath = os.path.abspath(os.curdir)
prefix = "hbp-00009"

folders = os.listdir(basepath)
folders = [el for el in folders if os.path.isdir(el) and not el.startswith(".")]
folders.remove("AP_propagation")
ap_path = os.path.join(basepath, "AP_propagation")
ap_file_list = os.listdir(ap_path)
print("working on {} folders".format(len(folders)))
for folder in folders:
    channel = folder.replace(".", "")
    folder_path = os.path.join(basepath, folder)
    folder_ext = "_".join([prefix, folder])
    folder_path_new = os.path.join(basepath, folder_ext)
    subfolders = os.listdir(folder_path)
    subfolders = [el for el in subfolders if os.path.isdir(os.path.join(folder_path,el)) and not el.startswith(".")]
    print("working on {} subfolders in\n{}\nrenamed to".format(len(subfolders), folder_path, folder_path_new))
    for subfolder in subfolders:
        subfolder_path = os.path.join(folder_path, subfolder)
        subfolder_ext = "__".join([folder_ext, subfolder])
        subfolder_path_new = os.path.join(folder_path, subfolder_ext)
        print("working in subfolder {}\nrenaming to {}".format(subfolder_path, subfolder_path_new))
        files = os.listdir(subfolder_path)
        files = [el for el in files if not el.startswith(".")]
        # We rename the files
        for file in files:
            fname, fext = os.path.splitext(file)
            name_parts = [subfolder_ext]
            name_parts.extend([el for el in fname.split("_") if not el in subfolder_ext])
            new_name = "_".join(name_parts)
            new_name = new_name + fext

            print("old_name", os.path.join(subfolder_path, fname))
            print("new_name", new_name)
            os.rename(os.path.join(subfolder_path, fname + fext), os.path.join(subfolder_path, new_name))

        # Now we rename the subfolder
        os.rename(subfolder_path, subfolder_path_new)

    # Identify corresponding files in the AP_propagation folder
    ap_files = filter(lambda el: el.startswith(channel), ap_file_list)
    if channel == "Kv101":
        ap_files = []
    if ap_files:
        ap_folder_ext = "_".join([folder_ext, "AP_propagation"])
        ap_path_new = os.path.join(folder_path, ap_folder_ext)
        os.mkdir(ap_path_new)
    for file in ap_files:
        old = os.path.join(ap_path, file)
        new = os.path.join(ap_path_new, ap_folder_ext + file.replace(channel,""))
        subprocess.call(["mv", old, new])

    # Now we rename the folder
    os.rename(folder_path, folder_path_new)
