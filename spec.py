import tables
import os
import utilities.misc
import numpy
import scipy.stats
import enum
import copy
import spc

class operation(enum.Enum):
    add = 1
    subtr = 2
    mult = 3
    div = 4

class Spec:
    """
        Spectrum

        id,
        x,
        data

    """
    #id of spec, used when loading from specs with no internal id...
    id = 0

    #create empy spec
    def __init__(self, id=0):
        if (id != 0):
            self.id = id
        else:
            Spec.id -= 1
            self.id = Spec.id

        self.wavelengths = []
        self.data = []

    def rebin(self, start, end, step=-1):
        # hopefully we have a constant step...
        if step == -1:
            step = numpy.average(numpy.diff(self.wavelengths))
        self.wavelengths, self.data = utilities.misc.rebin(self.wavelengths, self.data, start, end, step)

    def reverse(self):
        self.wavelengths = self.wavelengths[::-1]
        self.data = self.data[::-1]

    @staticmethod
    def loadSpecFromASCII(filename, delim=' '):
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
        xvalues = []
        yvalues = []
        for line in lines:
            split = line.split(delim)
            xvalues.append(float(split[0]))
            yvalues.append(float(split[1]))
        spc = Spec()
        spc.wavelengths = xvalues
        spc.data = yvalues
        return spc

    @staticmethod
    def loadSpecFromHDF5(filename):
        f = tables.open_file(filename)
        print(f.root.spec.id[0])
        spc = Spec(f.root.spec.id[0])
        spc.wavelengths = f.root.spec.q[:]
        spc.data = f.root.spec.data[:]
        f.close()
        return spc

    @staticmethod
    def loadAllSpecFromHDF5(path):
        spc = []
        for root, dirs, files in os.walk(path):
            for file in files:
                if file.endswith(".h5"):
                    spc.append(Spec.loadSpecFromHDF5(path + file))
        return spc

    @staticmethod
    def loadSpecFromSPC(filename):
        s = spc.File(filename)
        lines = s.data_txt().split("\n")

        #
        # FIXME: BAD BAD BAD BAD
        #
        id = -1
        try:
            fname = os.path.basename(filename)
            fname = fname.split('.')
            fname = fname[0].split('_')
            id = int(fname[-1])
        except Exception as err:
            print("Failed to parse filename for ID:", err)
            pass

        if id != -1:
            s = Spec(id)
        else:
            s = Spec()

        for l in lines[1:]:
            l = l.rstrip('\t\r\n')
            l = l.lstrip('\t\r\n')
            if l:
                x, y = l.split("\t")
                s.wavelengths.append(float(x))
                s.data.append(float(y))

        return s

    @staticmethod
    def loadAllSpecFromSPC(path):
        spc = []
        for root, dirs, files in os.walk(path):
            for file in files:
                if file.endswith(".spc"):
                    spc.append(Spec.loadSpecFromSPC(path + file))
        return spc

def getSpecByID(specs, id):
    s = []
    for spec in specs:
        if isinstance(id, list):
            if spec.id in id:
                s.append(spec)
        else:
            if spec.id == id:
                s.append(spec)
    if len(s) == 1:
        return s[0]
    else:
        return s

def specOperation(inputSpecs, secondOperand, op):
    if isinstance(inputSpecs, list):
        specs = []
        for s in inputSpecs:
            specs.append(specOperation(s, secondOperand, op))
        return specs

    specs = copy.deepcopy(inputSpecs)

    if isinstance(secondOperand, Spec):
        overlap = False
        overlapMin = numpy.maximum(numpy.amin(specs.wavelengths), numpy.amin(secondOperand.wavelengths))
        overlapMax = numpy.minimum(numpy.amax(specs.wavelengths), numpy.amax(secondOperand.wavelengths))
        if overlapMax > overlapMin:
            overlap = True
        if not overlap:
            raise Exception("Specs do not overlap.")

        step = numpy.average(numpy.diff(specs.wavelengths))

        specs.rebin(overlapMin, overlapMax, step)
        secondOperand.rebin(overlapMin, overlapMax, step)

        if op == operation.add:
            specs.data = numpy.add(specs.data, secondOperand.data)
        elif op == operation.subtr:
            specs.data = numpy.subtract(specs.data, secondOperand.data)
        elif op == operation.mult:
            specs.data = numpy.multiply(specs.data, secondOperand.data)
        elif op == operation.div:
            specs.data = numpy.divide(specs.data, secondOperand.data)

    else:
        if op == operation.add:
            specs.data = numpy.add(specs.data, secondOperand)
        elif op == operation.subtr:
            specs.data = numpy.subtract(specs.data, secondOperand)
        elif op == operation.mult:
            specs.data = numpy.multiply(specs.data, secondOperand)
        elif op == operation.div:
            specs.data = numpy.divide(specs.data, secondOperand)

    return specs


def getSpecInstrumentationAbsError(spec, minDark=220, maxDark=300):
    '''Absolute error from instrumentation noise'''
    darkNoise = Spec(id=-1)
    darkNoise.data = spec.data
    darkNoise.wavelengths = spec.wavelengths
    darkNoise.rebin(minDark, maxDark)
    return numpy.std(darkNoise.data)

def getSpecInstrumentationRelError(spec, minDark=220, maxDark=300):
    '''Relative error from instrumentation noise'''
    absErr = getSpecInstrumentationAbsError(spec, minDark, maxDark)
    return numpy.abs([absErr/i for i in spec.data])

def getSpecStudnetTError(specs, p=0.001):
    '''Absolute error from averaging'''

    avgSpec = numpy.average([spec.data for spec in specs], axis = 0)
    #Standard Error  = Std. Deviation / sqrt(N)
    specStandardError = numpy.divide(numpy.std([numpy.subtract(spec.data, avgSpec) for spec in specs], axis=0), numpy.sqrt(len(specs)))

    tError = numpy.multiply(specStandardError,  scipy.stats.t.ppf(p, len(specs) - 1))

    return numpy.abs(tError)

def getSpecStudentReTError(specs, p=0.001):
    '''Relative error from averaging'''
    tError = getSpecStudnetTError(specs, p)
    avgSpec = numpy.average([spec.data for spec in specs], axis = 0)
    return numpy.divide(tError, avgSpec)

def getSpecAbsError(spec, darkMin=220, darkMax=300, studentErr=None):
    specInstrumentationErr = getSpecInstrumentationAbsError(spec, darkMin, darkMax)
    if studentErr is None:
        return [specInstrumentationErr for i in spec.data]
    else:
        return numpy.sqrt(numpy.add(numpy.multiply(studentErr, studentErr), numpy.multiply(specInstrumentationErr, specInstrumentationErr)))

def getSpecRelError(theSpec, darkMin=220, darkMax=300, studentErr=None):
    err = getSpecAbsError(theSpec, darkMin=darkMin, darkMax=darkMax, studentErr=studentErr)
    return numpy.abs(numpy.divide(err, theSpec.data))