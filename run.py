import os

dataset = "deep"
vecdim  = 96
size_base = 1
format = "Float"

path_dataset = "/home/usr-xkIJigVq/DataSet/" + dataset
path_datasize = path_dataset + "/" + dataset + str(size_base) + "m"

path_base   = path_datasize + "/base." + str(size_base) + "m.fbin"
path_gt     = path_datasize + "/groundtruth." + str(size_base) + "m.bin"
path_query  = path_dataset + "/query.public.10K.fbin"

path_index = "graphindex/" + dataset + str(size_base) + "m"

def main():
    cmd_make = "cd build && cmake .. && make -j"
    os.system(cmd_make)

    command = "./Release/indexbuilder " + \
                "-c buildconfig.ini " + \
                "-d " + str(vecdim) + " " + \
                "-v " + format + " " + \
                "-f DEFAULT " + \
                "-i " + path_base + " " + \
                "-o " + path_index + " " + \
                "-a SPANN"
    os.system(command)

main()

# "-c buildconfig.ini " + \