import numpy
import itertools
import scipy.interpolate
import scipy.signal

__author__ = 'Kryosugarra'

def rebin(xValues, yValues, start, end, step):
    assert start >= xValues[0]
    assert end <= xValues[-1]
    assert step != 0
    f = scipy.interpolate.interp1d(xValues, yValues)
    newX = numpy.linspace(start, end, (end-start+step)/(step))
    newY = []
    for x in newX:
        newY.append(f(x))
    return newX, newY


def rebinFile(filename, newfile, start, end, step):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    xvalues = []
    yvalues = []
    for line in lines:
        split = line.split()
        xvalues.append(float(split[0]))
        yvalues.append(float(split[1]))
    xvalues, yvalues = rebin(xvalues, yvalues, start, end, step)
    f = open(newfile, 'w')
    for i in range(len(xvalues)):
        f.write('{0:15f} \t {1:15f} \n'.format(float(xvalues[i]), float(yvalues[i])))
    f.close()