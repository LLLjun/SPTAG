import os

dataset = "deep"
vecdim  = 96
size_base = 1
format = "Float"

stage = "search"
cmake_mode = "Release"
alg = "SPANN"

path_dataset = "/home/usr-xkIJigVq/DataSet/" + dataset
path_datasize = path_dataset + "/" + dataset + str(size_base) + "m"

path_base   = path_datasize + "/base." + str(size_base) + "m.fbin"
path_gt     = path_datasize + "/groundtruth." + str(size_base) + "m.bin"
path_query  = path_dataset + "/query.public.10K.fbin"

path_index = "graphindex/" + alg + "/" + dataset + str(size_base) + "m_test"

def main():
    cmd_make = "cd build && cmake -DCMAKE_BUILD_TYPE=" + cmake_mode + " .."
    if stage == "both" or stage == "build":
        cmd_make += " && make indexbuilder -j"
    if stage == "both" or stage == "search":
        cmd_make += " && make indexsearcher -j"
    os.system(cmd_make)

    cmd_build = "./Release/indexbuilder " + \
                "-c buildconfig.ini " + \
                "-d " + str(vecdim) + " " + \
                "-v " + format + " " + \
                "-f DEFAULT " + \
                "-i " + path_base + " " + \
                "-o " + path_index + " " + \
                "-a " + alg + " " + \
                "-t 112"
    if stage == "both" or stage == "build":
        os.system(cmd_build)

    for t in [8]:
        cmd_search = "./Release/indexsearcher " + \
                    "-d " + str(vecdim) + " " + \
                    "-v " + format + " " + \
                    "-i " + path_query + " " + \
                    "-r " + path_gt + " " + \
                    "-f DEFAULT " + \
                    "-t " + str(t) + " " + \
                    "-k 10 " + \
                    "-x " + path_index + " " + \
                    "-m 1#4#8#16#32#1024#2048#4096#8192"
        if stage == "both" or stage == "search":
            os.system(cmd_search)

main()

# "-c buildconfig.ini " + \