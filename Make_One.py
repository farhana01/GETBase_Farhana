import os
tld = "./Output/"


def main():
    best_idxs = ["01", "06", "11", "16", "21", "26"]

    folds = os.listdir(tld)
    for fold in folds:
        cnt = 0
        with open(tld + fold + "/best.dat", "w") as out:
            for idx in best_idxs:
                with open(tld + fold + "/best" + idx + ".dat", "r") as inn:
                    lines = inn.readlines()
                    for line in lines:
                        out.writelines(line)
                    pass
                pass
            pass
        pass

    pass


main()