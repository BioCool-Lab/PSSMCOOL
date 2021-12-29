import pathlib
import os

current_path = str(pathlib.Path().resolve())
run_time_folder = current_path + '/runtime_files'
files_list = os.listdir(run_time_folder)
final_list = []
for name in files_list:
    txt_file_path = run_time_folder + '/' + name
    with open(txt_file_path)as txt_file:
        for line in txt_file:
            if "The specified output file name is" in line:
                protein_name = line.strip().split('rpssm_')[1]
            if "command took" in line:
                run_time = line.split('(')[1].split(' total')[0]
        final_list.append([protein_name, run_time])

with open(current_path+'/runtimes.csv', 'w')as csv_write_file:
    csv_write_file.write('' + ',' + "Proteins' name" + ',' + 'Runtime' + '\n' )
    counter = 0
    for i in final_list:
        counter += 1
        csv_write_file.write(str(counter) + ',' + i[0] + ',' + i[1] + '\n')
