def get_data(pathToFile):
    phenotypeFreq = {}
    with open(pathToFile) as f:
        while (True):
            line = f.readline()

            if not line:
                break

            data = line.split(" ")
            phenotypeFreq[data[0]] = int(data[1].replace("\n", ""))
    return phenotypeFreq
