#!/usr/bin/env python

import unittest
from gamefft.test.tests import ResultValidationTest

import sys,os
sys.path.remove(os.getcwd())


if __name__=="__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(ResultValidationTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
