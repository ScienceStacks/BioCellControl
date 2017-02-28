import unittest
from chemotaxis_model import ChemotaxisModel, ReceptorStates, State


IGNORE_TEST = False


#############################
# Tests
#############################
# pylint: disable=W0212,C0111,R0904
class TestChemotaxisModelBasic(unittest.TestCase):

  def setUp(self):
    self.model = ChemotaxisModel()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue('Yp =' in self.model._model)

  def testInitialize(self):
    if IGNORE_TEST:
      return
    rr = self.model.initialize()
    self.assertEqual(rr.R, 0.3e-6)

  def testRun(self):
    if IGNORE_TEST:
      return
    self.model.initialize()
    self.model.run()
    self.assertTrue(abs(self.model._rr.R - 0.3e-6) <= 1e-8)


# pylint: disable=W0212,C0111,R0904
class TestState(unittest.TestCase):

  def setUp(self):
    self.model = ChemotaxisModel()
    self.model.initialize()
    self.model.run(start=0, end=10)

  def _testBasics(self, l, p, m, expected_name):
    """
    :param bool l: is ligand bound
    :param bool p: is phosphorylated
    :param int m: methylation 
    """
    state = State(l, p, m, self.model.getResult())
    name = state.getName()
    self.assertEqual(name, expected_name)
    pairs = zip(state.getData(), self.model.getResult()[name])
    b = all([abs(x-y) <= 1e-3*abs(min(x,y)) for x,y in pairs])
    if not b:
      import pdb; pdb.set_trace()
    self.assertTrue(b)

  def testBasics(self):
    self._testBasics(True, True, 2, '[LT2p]')
    self._testBasics(False, True, 2, '[T2p]')
    self._testBasics(True, False, 2, '[LT2]')
    self._testBasics(True, True, 3, '[LT3p]')


# pylint: disable=W0212,C0111,R0904
class TestReceptorStates(unittest.TestCase):

  def setUp(self):
    self.model = ChemotaxisModel()
    self.model.initialize()
    self.model.run(start=0, end=10)

if __name__ == '__main__':
  unittest.main()
