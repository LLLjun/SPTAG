import os

class RunConfig:
    dataset = None
    vecdim  = None
    format  = None
    qsize   = None
    path_gt     = None
    path_base   = None
    path_query  = None
    path_index  = None
    path_result = None
    path_config = "Config/billion/"

    def __init__(self, dataname, mode="billion"):
        self.dataset = dataname
        self.mode = mode
        self.CheckDataset()

    def CheckDataset(self):
        if self.dataset == "sift":
            self.vecdim  = 128
            self.format = "UInt8"
            self.qsize   = 10000
        elif self.dataset == "spacev":
            self.vecdim  = 100
            self.format = "Int8"
            self.qsize   = 29316
        elif self.dataset == "deep":
            self.vecdim  = 96
            self.format = "Float"
        elif self.dataset == "turing":
            self.vecdim  = 100
            self.format = "Float"
        else:
            print("Error, unexisted dataname")
            exit(1)

        if self.mode == "billion":
            path_dataset = "../dataset/billion/" + self.dataset
            self.path_gt     = path_dataset + "/groundtruth.bin"
            self.path_base   = path_dataset + "/base"
            self.path_query  = path_dataset + "/query"
            self.path_index  = "graphindex/" + self.dataset
        elif self.mode == "toy":
            path_dataset = "../dataset/sift"
            self.path_gt     = path_dataset + "/sift10m/groundtruth.10m.bin"
            self.path_base   = path_dataset + "/sift10m/base.10m.u8bin"
            self.path_query  = path_dataset + "/query.public.10K.u8bin"
            self.path_index  = "graphindex/toy/sift10m_12"
            self.path_result = "output/result/toy/sift10m_rc10.log"

    def make(self, stage):
        cmake_mode = "Release"
        cmd_make = "cd build && cmake -DCMAKE_BUILD_TYPE=" + cmake_mode + " .."
        if stage == "build":
            cmd_make += " && make indexbuilder -j"
        elif stage == "search":
            cmd_make += " && make indexsearcher -j"
        os.system(cmd_make)

    def run(self, stage):
        alg = "SPANN"
        config_file = self.path_config + self.dataset + ".ini"
        if stage == "build":
            cmd_build = "./Release/indexbuilder " + \
                        "-c " + " " + config_file + \
                        "-d " + str(self.vecdim) + " " + \
                        "-v " + self.format + " " + \
                        "-f DEFAULT " + \
                        "-i " + self.path_base + " " + \
                        "-o " + self.path_index + " " + \
                        "-a " + alg + " " + \
                        "-t 200"
            os.system(cmd_build)
        elif stage == "search":
            batch_size = 100
            for t in [1]:
                cmd_search = "./Release/indexsearcher " + \
                            "-d " + str(self.vecdim) + " " + \
                            "-v " + self.format + " " + \
                            "-i " + self.path_query + " " + \
                            "-r " + self.path_gt + " " + \
                            "-f DEFAULT " + \
                            "-t " + str(t) + " " + \
                            "-k 10 " + \
                            "-x " + self.path_index + " " + \
                            "-m 100 " + \
                            "-b " + str(batch_size)
                            # "-o " + self.path_result
                os.system(cmd_search)

    def make_run(self, stage):
        self.make(stage)
        self.run(stage)


def main():
    dataname = "sift"
    stage = "search"
    mode = "billion"

    Run = RunConfig(dataname, mode)
    Run.make_run(stage)

if __name__ == '__main__':
    main()