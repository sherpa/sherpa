#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2013)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
from sherpa.utils import *
from sherpa.utils.err import SherpaErr
from sherpa.utils.err import ModelErr

class test_err(SherpaTestCase):

   def test_NewSherpaErr(self):
	class OldSherpaErr( Exception ): 
     		"Old class for all Sherpa exceptions" 
     		def __init__( self, dict, key, *args ): 
         		if dict.has_key( key ): 
             			errmsg = dict[ key ] % args 
         		else: 
             			errmsg = "unknown key '%s'" % key 
        	  	Exception.__init__(self, errmsg) 
        
	dict = {'simple':'simple message', 'arg':'argument: %s'}
	
	# Test 1: verify that a correct call of the new constructor has the same result of the old one
	err = SherpaErr(dict, 'simple')
	old_err = OldSherpaErr(dict, 'simple')
        self.assertEqual(err.message, old_err.message)

	# Test 2: same as before, but with string placeholders
	err = SherpaErr(dict, 'arg', 'foo')
	self.assertEqual('argument: foo', err.message)

	# Test #3: verify that a call without a key results in a generic message being produced
	err = SherpaErr(dict)
	self.assertEqual('Generic Error', err.message)
	
	# Test #4: verify the user's expected behavior, i.e. a string is provided as error message
	err = SherpaErr(dict, 'My Error')
	self.assertEqual('My Error', err.message)
	
	# Test #5: verify the user provided example, which exercises a derived class
	err = ModelErr("Unable to frobnicate model %s" % 'modelname') 
        self.assertEqual('Unable to frobnicate model modelname', err.message)


if __name__ == '__main__':

    import sherpa.utils as utils
    SherpaTest(utils).test()
