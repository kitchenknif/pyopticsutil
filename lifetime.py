import os
import numpy
import scipy.optimize

from .utilities.geometry import get1DPeak


def load_LT_from_path(path):
    prefix = "G2data"
    folderprefix = "LTdata"
    # param1 = "wire"
    # param2 = "pos"
    # metadatapostfix = "inf"
    extension = ".dat"
    filepairs = []

    for root, dirs, files in os.walk(os.path.normpath(path)):
        if folderprefix in root:
            for file in files:
                if prefix in file and extension in file:
                    if not "inf" in file:
                        parsed = file.split(".")
                        (header,) = [s for s in files if parsed[0] + "inf." + parsed[1] in s]
                        filepairs.append([os.path.join(root, file), os.path.join(root, header)])
    LTs = []
    for file, header in filepairs:
        LT = read_LT(file)
        com = read_header(header)
        LTs.append([LT, com])
    return LTs


def read_header(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    comment = ""
    for line in lines:
        if line.startswith("Comment"):
            comment = line[8:]
    return comment


def read_LT(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    time = []
    intensity = []
    for line in lines:
        s = line.split('\t')
        time.append(float(s[1]))
        intensity.append(float(s[3]))
    return time, intensity


def filter_LT(t, i, *, binsize=5):
    x, y = get1DPeak(t, i, interpolate=False, indicies=True, use="mode")

    #Y
    base = numpy.average(i[:x])
    i = numpy.subtract(i, base)
    i = numpy.divide(i[x:], numpy.max(i[x:]))
    #i = numpy.log(i)

    #X
    t = t[x:]
    t = numpy.subtract(t, t[0])

    i_f = []
    i_f_rms = []
    t_f = []
    for j in range(0, len(i), binsize):
        bin = i[j: j+binsize]
        if j > 0 and (numpy.average(bin) < 0 or numpy.average(bin) <= numpy.std(bin)):
            break
        i_f.append(numpy.average(bin))
        i_f_rms.append(numpy.std(bin))
        bin_t = t[j: j+binsize]
        t_f.append(numpy.average(bin_t))

    return t_f, i_f, i_f_rms


def fit_LT(t, i):
    def piecewise_linear(x, x0, y0, k1, k2):
        return numpy.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

    params, pcov = scipy.optimize.curve_fit(piecewise_linear, t, numpy.log(i))
    err = numpy.sqrt(numpy.diag(pcov))

    return params[-2], params[-1], err, piecewise_linear(t, *params)