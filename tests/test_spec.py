from unittest import TestCase
import spec
__author__ = 'Kryosugarra'


class TestSpec(TestCase):
    def test_rebin(self):
        self.fail()

    def test_loadSpecFromASCII(self):
        self.fail()

    def test_loadSpecFromHDF5(self):
        self.fail()

    def test_loadAllSpecFromHDF5(self):
        self.fail()

    def test_loadSpecFromSPC(self):
        s = spec.Spec.loadSpecFromSPC("Data_17_03_51791.spc")
        print(s.id, len(s.data), len(s.wavelengths))

        if s.id != 51791:
            self.fail()

    def test_loadAllSpecFromSPC(self):
        self.fail()
