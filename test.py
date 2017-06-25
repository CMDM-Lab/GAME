#!/usr/bin/env python
import sys,os
sys.path.remove(os.getcwd())


import unittest
from gamefft.test.tests import ResultValidationTest


if __name__=="__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(ResultValidationTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
