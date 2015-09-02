# -*- coding: utf-8 -*-

import numpy
import itertools
import scipy.interpolate
import scipy.signal

__author__ = 'Kryosugarra'


#
# Peak find
#
def get1DPeak(x, y, *, inverse=False, interpolate=True, indicies=False, use='mass'):
    assert len(x) == len(y)
    if inverse:
        maxY = numpy.max(y)
        y = numpy.multiply(y, -1)
        y = numpy.add(y, maxY)
    peakPos = -1
    if use == 'mass':
        peakPos = numpy.divide(numpy.sum(numpy.multiply(x, y)), numpy.sum(y))
    elif use == 'mode':
        peakPos = x[numpy.where(y == numpy.max(y))[0][0]]
    elif use == 'median':
        mode = numpy.median(y)
        diffs = [numpy.abs(i - mode) for i in y]
        peakPos = x[numpy.where(diffs == numpy.min(diffs))[0][0]]
    elif use == 'average':
        avg = numpy.average(y)
        diffs = [numpy.abs(i - avg) for i in y]
        peakPos = x[numpy.where(diffs == numpy.min(diffs))[0][0]]
    else:
        raise Exception('Bad peak find mechanism')

    if indicies:
        diffs = [i - peakPos for i in x]
        peakIndex = diffs.index(numpy.min(numpy.abs(diffs)))
        peakPos = peakIndex
        peakVal = y[peakIndex]
    elif interpolate:
        f = scipy.interpolate.interp1d(x, y)
        peakVal = f(peakPos)
    else:
        diffs = [i - peakPos for i in x]
        peakIndex = diffs.index(numpy.min(numpy.abs(diffs)))
        peakPos = x[peakIndex]
        peakVal = y[peakIndex]

    return peakPos, peakVal


def get2DPeak(x, y, z, inverse=False, interpolate=True, indicies=False):
    if inverse:
        maxZ = numpy.max(z)
        z = numpy.multiply(z, -1)
        z = numpy.add(z, maxZ)

    peakX = -1
    peakY = -1
    peakVal = -1
    #if len(z) == len(x)*len(y):
    #    peakX = numpy.divide(numpy.sum(numpy.multiply(x, z)), numpy.sum(z))
    #    peakY = numpy.divide(numpy.sum(numpy.multiply(y, z)), numpy.sum(z))
    #el
    if z.shape == (len(x), len(y)):
        peakX = 0
        peakY = 0
        for i in range(len(x)):
            for j in range(len(y)):
                peakX += z[i, j]*x[i]
                peakY += z[i, j]*y[j]
        peakX = peakX/numpy.sum(z)
        peakY = peakY/numpy.sum(z)
        if indicies:
            diffs = [numpy.abs(i - peakX) for i in x]
            xIndex = diffs.index(numpy.min(diffs))
            peakX = xIndex
            diffs = [numpy.abs(i - peakY) for i in y]
            yIndex = diffs.index(numpy.min(diffs))
            peakY = yIndex
            peakVal = z[xIndex, yIndex]
        elif interpolate:
            f = scipy.interpolate.interp2d(x, y, z)
            peakVal = f(peakX, peakY)
        else:
            diffs = [numpy.abs(i - peakX) for i in x]
            xIndex = diffs.index(numpy.min(diffs))
            peakX = x[xIndex]
            diffs = [numpy.abs(i - peakY) for i in y]
            yIndex = diffs.index(numpy.min(diffs))
            peakY = y[yIndex]
            peakVal = z[xIndex, yIndex]
    else:
        raise Exception('Bad shape of x={0}, y={1}, z={2} params'.format(len(x), len(y), z.shape))

    return peakX, peakY, peakVal

def getPeakCurvatureRadius(x, y, peakMin = -1, peakMax = -1):
    print('FIXME')
    f = scipy.interpolate.interp1d(x, y, 'cubic')
    m = 0
    peakX = -1
    for x_i, y_i in itertools.zip_longest(x, y):
        if y_i >= m:
            m = y_i
            peakX = x_i

    #peakX = get1DPeak(x, y)

    #deriv1
    dx = numpy.diff(x)
    y1 = numpy.divide(numpy.diff(f(x)), dx)
    x1 = numpy.add(numpy.divide(numpy.diff(x), 2), x[0:-1])
    f1 = scipy.interpolate.interp1d(x1, y1, 'cubic')

    #deriv2
    dx = numpy.diff(x1)
    y2 = numpy.divide(numpy.diff(f1(x1)), dx)
    x2 = numpy.add(numpy.divide(numpy.diff(x1), 2), x1[0:-1])
    f2 = scipy.interpolate.interp1d(x2, y2, 'cubic')

    R = numpy.abs(numpy.divide(numpy.power(numpy.add(1, numpy.power(f1(x2), 2)), 3.0/2.0), f2(x2)))
    Rint = scipy.interpolate.interp1d(x2, R, 'cubic')

    #plt.plot(x, f(x))
    #plt.plot(x1, f1(x1))
    #plt.plot(x2, f2(x2))
    #plt.plot(x2, Rint(x2))

    numR = -1
    if peakMax != -1 and peakMin != -1:
        tehR, n = 0, 0
        for x, y in itertools.zip_longest(x2, R):
            if 1100 <= x <= 1125:
            #if 450 <= x <= 500:
                n += 1
                tehR += y
        numR = tehR/n

    return Rint(peakX), peakX, numR


#
# Array slicing
#
def getArraySlice(array, x0, y0, x1, y1):
    # !! x0, y0, x1, y1  are in _pixel_ coordinates !!

    length = int(numpy.hypot(x1-x0, y1-y0))
    x, y = numpy.linspace(x0, x1, length), numpy.linspace(y0, y1, length)

    # Extract the values along the line
    slice = array[x.astype(numpy.int), y.astype(numpy.int)]
    pixelSize = len(slice)
    return slice


def getSliceByAngle(z, angle, x0, y0, x=-1, y=-1, degrees=True):
    """
    returns slice of rectangular array passing through point (x0, y0) at an angle to a vertical line (x0, :)
    all coordinates in PIXELS
    """

    if degrees:
        angle = angle*numpy.pi/180

    while angle >= 2*numpy.pi:
        angle -= 2*numpy.pi
    while angle <= -2*numpy.pi:
        angle += 2*numpy.pi

    reverse = False
    if angle > numpy.pi/2:
        angle -= numpy.pi
        reverse = True
    if angle < -numpy.pi/2:
        angle += numpy.pi
        reverse = True


    if x == -1:
        x = range(len(z[:][0]))
    if y == -1:
        y = range(len(z[0][:]))


    #start
    x1, y1 = x0, 0
    x2, y2 = x0, len(y) - 1
    if numpy.abs(angle) < numpy.arctan(1/numpy.max([len(x), len(y)])):
        x1, y1 = x0, 0
        x2, y2 = x0, len(y) - 1
    elif numpy.abs(numpy.abs(angle) - numpy.pi/2) >= numpy.abs(numpy.arctan(1/numpy.max([len(x), len(y)]))):
        x1 = x0 + y0 * numpy.tan(angle)
        x2 = x0 - (len(y) - y0 - 1) * numpy.tan(angle)
        #X out-of-bounds
        if x1 >= len(x):
            y1 = (x1 - (len(x) - 1))/numpy.tan(numpy.abs(angle))
            x1 = len(x) - 1
        elif x1 < 0:
             y1 = numpy.abs(x1)/numpy.tan(numpy.abs(angle))
             x1 = 0

        if x2 < 0:
            y2 = len(y) - 1 - numpy.abs(x2)/numpy.tan(numpy.abs(angle))
            x2 = 0
        elif x2 >= len(x):
            y2 = len(y) - 1 - (x2 - (len(x) - 1))/numpy.tan(numpy.abs(angle))
            x2 = len(x) - 1

    else:
        x1, y1 = 0, y0
        x2, y2 = len(x) - 1, y0
    if not reverse:
        return getArraySlice(z, x1, y1, x2, y2), x1, y1, x2, y2
    else:
        return getArraySlice(z, x2, y2, x1, y1), x2, y2, x1, y1

