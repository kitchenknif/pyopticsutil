import tables
import os
import utilities.misc
import numpy
import scipy.stats
import enum
import copy
import spc.spc as spc

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

        self.x = numpy.empty(1, dtype=numpy.float64)
        self.wavelengths = self.x  # Time to deprecate wavelengths, since use is wider.
        self.data = numpy.empty(1, dtype=numpy.float64)

    def rebin(self, start, end, step=-1):
        # hopefully we have a constant step...
        if step == -1:
            step = numpy.average(numpy.diff(self.wavelengths))
        self.wavelengths, self.data = utilities.misc.rebin(self.wavelengths, self.data, start, end, step)

    def reverse(self):
        self.wavelengths = self.wavelengths[::-1]
        self.data = self.data[::-1]

    def save_as_ascii(self, filename):
        with open(filename, 'w') as f:
            f.write("#X (nm), Y\n")
            for x, y in zip(self.wavelengths, self.data):
                f.write("{}, {}\n".format(x, y))

    @staticmethod
    def specFromXYData(x, y):
        spc = Spec()
        spc.wavelengths = numpy.asarray(x)
        spc.data = numpy.asarray(y)
        return spc

    @staticmethod
    def loadSpecFromASCII(filename, *, decimal='.', delim=' ', headerlines=0):
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
        xvalues = []
        yvalues = []
        for line in lines[headerlines:]:
            if not line.startswith('#'):
                line = line.replace(decimal, '.')
                split = line.split(delim)
                xvalues.append(float(split[0]))
                yvalues.append(float(split[1]))
        spc = Spec()
        spc.wavelengths = numpy.asarray(xvalues)
        spc.data = numpy.asarray(yvalues)
        spc.filename = os.path.basename(filename)
        return spc

    @staticmethod
    def loadAllSpecFromASCII(path, *, extension='.txt', decimal='.', delim=' ', headerlines=0):
        spc = []
        for root, dirs, files in os.walk(path):
            for file in files:
                if file.endswith(extension):
                    spc.append(Spec.loadSpecFromASCII(os.path.join(root, file),
                                                      decimal=decimal, delim=delim, headerlines=headerlines))
        return spc


    @staticmethod
    def loadSpecFromHDF5(filename):
        f = tables.open_file(filename)
        print(f.root.spec.id[0])
        spc = Spec(f.root.spec.id[0])
        spc.wavelengths = numpy.asarray(f.root.spec.q[:])
        spc.data = numpy.asarray(f.root.spec.data[:])
        spc.filename = os.path.basename(filename)
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

        xvals = []
        yvals = []
        for l in lines[1:]:
            l = l.rstrip('\t\r\n')
            l = l.lstrip('\t\r\n')
            if l:
                x, y = l.split("\t")
                xvals.append(float(x))
                yvals.append(float(y))

        s.wavelengths = numpy.asarray(xvals)
        s.data = numpy.asarray(yvals)
        s.filename = os.path.basename(filename)
        return s

    @staticmethod
    def loadAllSpecFromSPC(path):
        spc = []
        for root, dirs, files in os.walk(path):
            for file in files:
                if file.endswith(".spc"):
                    spc.append(Spec.loadSpecFromSPC(os.path.join(root, file)))
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

    sp = copy.deepcopy(inputSpecs)

    if isinstance(secondOperand, Spec):
        overlap = False
        overlapMin = numpy.maximum(numpy.amin(sp.wavelengths), numpy.amin(secondOperand.wavelengths))
        overlapMax = numpy.minimum(numpy.amax(sp.wavelengths), numpy.amax(secondOperand.wavelengths))
        if overlapMax > overlapMin:
            overlap = True
        if not overlap:
            raise Exception("Specs do not overlap.")

        step = numpy.average(numpy.diff(sp.wavelengths))

        sp.rebin(overlapMin, overlapMax, step)
        secondOp = copy.deepcopy(secondOperand)
        secondOp.rebin(overlapMin, overlapMax, step)

        if op == operation.add:
            sp.data = numpy.add(sp.data, secondOp.data)
        elif op == operation.subtr:
            sp.data = numpy.subtract(sp.data, secondOp.data)
        elif op == operation.mult:
            sp.data = numpy.multiply(sp.data, secondOp.data)
        elif op == operation.div:
            sp.data = numpy.divide(sp.data, secondOp.data)

    else:
        if op == operation.add:
            sp.data = numpy.add(sp.data, secondOperand)
        elif op == operation.subtr:
            sp.data = numpy.subtract(sp.data, secondOperand)
        elif op == operation.mult:
            sp.data = numpy.multiply(sp.data, secondOperand)
        elif op == operation.div:
            sp.data = numpy.divide(sp.data, secondOperand)

    return sp

def getSpecBaseline(sp, poly=1):
    s = Spec()
    coeffs = numpy.polyfit(sp.wavelengths, sp.data, poly)
    p = numpy.poly1d(coeffs)
    s.wavelengths = copy.copy(sp.wavelengths)
    s.data = numpy.asarray([p(lambd) for lambd in s.wavelengths])
    return s

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