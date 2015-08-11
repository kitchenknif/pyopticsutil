# -*- coding: utf-8 -*-

import numpy

def loadSpecFromTXT(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    xvalues = []
    yvalues = []
    for line in lines:
        split = line.split()
        xvalues.append(float(split[0]))
        yvalues.append(float(split[1]))
    return xvalues, yvalues

def loadIAPro2DArray(filename):
    f = open(filename)
    lines = f.readlines()
    header = [lines[i] for i in range(0, 11)]
    header = [s.split('=') for s in header]
    xstep = 1000
    xsize = 1
    ystep = 1000
    ysize = 1
    zstep = 1000
    for s in header:
        if s[0] == "XSize":
            xstep *= float(s[1])
        elif s[0] == "YSize":
            ystep *= float(s[1])
        elif s[0] == "XPoints":
            xsize = int(s[1])
            xstep /= float(s[1])
        elif s[0] == "YPoints":
            ysize = int(s[1])
            ystep /= float(s[1])
        elif s[0] == "ZUnits":
            if 'nm' in s[1]:
                zstep = 1
            elif 'µm' in s[1]:
                zstep = 1000
        elif s[0] == "XUnits":
            if 'nm' in s[1]:
                xstep /= 1000
            elif 'µm' in s[1]:
                xstep /= 1
        elif s[0] == "YUnits":
            if 'nm' in s[1]:
                ystep /= 1000
            elif 'µm' in s[1]:
                ystep /= 1
    data = [s.split() for s in lines[11:]]
    x = [i*xstep for i in range(xsize)]
    y = [i*ystep for i in range(ysize)]
    z = numpy.zeros((xsize, ysize))
    for p in data:
        z[int(p[0])][int(p[1])] = (float(p[2])*zstep)

    return x, xstep, y, ystep, z, zstep


def loadIAPro1DArray(filename, xstep=1000, ystep=1000):
    f = open(filename)
    lines = f.readlines()
    header = lines[0]
    data = [s.split() for s in lines[1:]]
    x = []
    y = []
    for p in data:
        x.append(float(p[0])*xstep)
        y.append(float(p[1])*ystep)
    return x, y