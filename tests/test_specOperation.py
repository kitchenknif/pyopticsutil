from unittest import TestCase
import spec

__author__ = 'Kryosugarra'


class TestSpecOperation(TestCase):
    def test_specOperation(self):
        s = spec.Spec()
        s.data = [1, 2, 3, 4]
        s.wavelengths = [1, 2, 3, 4]

        s2 = spec.Spec()
        s2.data = [1, 2, 3, 4]
        s2.wavelengths = [1, 2, 3, 4]

        res = spec.specOperation(s, s2, spec.operation.subtr)

        if res == s:
            self.fail("Result should be copy, not same spec")

        if res.id != s.id:
            self.fail("Result id not equal to first operand")

        for i in res.data:
            if i != 0:
                self.fail("Math failed")

        for i in range(len(res.wavelengths)):
            if res.wavelengths[i] != s.wavelengths[i]:
                self.fail("rebin failed")

    def test_specOperation_List(self):
        self.fail()