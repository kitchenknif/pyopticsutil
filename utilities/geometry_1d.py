# -*- coding: utf-8 -*-

import numpy
import itertools
import scipy.interpolate
import scipy.signal

import numpy
import itertools
import scipy.interpolate
import scipy.signal
import scipy.optimize

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
    elif use == 'max':
        max = numpy.max(y)
        diffs = [numpy.abs(i - max) for i in y]
        peakPos = x[numpy.where(diffs == numpy.min(diffs))[0][0]]
    elif use == 'FWHW':
        max = numpy.max(y)
        #center = y.index(max)
        center = numpy.where(y==max)[0][0]
        diffs = [numpy.abs(i - max/2.) for i in y]

        diffs_r = diffs[center:]
        diffs_l = diffs[:center]
        firstIntersection = diffs_l.index(numpy.min(numpy.abs(diffs_l)))

        secondIntersectoin = center + diffs_r.index(numpy.min(numpy.abs(diffs_r)))

        peakPos = x[firstIntersection] /2. + x[secondIntersectoin]/2.
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
        diffs = [numpy.abs(i - peakPos) for i in x]
        #peakIndex = diffs.index(numpy.min(numpy.abs(diffs)))
        peakIndex = diffs.index(numpy.min(diffs))
        peakPos = x[peakIndex]
        peakVal = y[peakIndex]

    return peakPos, peakVal


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
# Function fitting
#
def gaussian(x,a,x0,sigma):
    return a*numpy.exp(-(x-x0)**2/(2*sigma**2))


def fit_gaussian(s, bounds=None, method="lm"):
    pos, int = get1DPeak(s.wavelengths, s.data, interpolate=True, indicies=False, use="mass")
    sigma = numpy.sqrt(numpy.sum(numpy.multiply(s.data, numpy.subtract(s.wavelengths, pos)**2))/len(s.data))

    if bounds:
        popt, pcov = scipy.optimize.curve_fit(gaussian, s.wavelengths, s.data,
                               p0=[int, pos, sigma], method=method,
                               bounds=bounds)
    else:
        popt, pcov = scipy.optimize.curve_fit(gaussian, s.wavelengths, s.data,
                               p0=[int, pos, sigma], method=method)
    return popt