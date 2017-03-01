import unittest
from chemotaxis_model import ChemotaxisModel, ReceptorStates, State,  \
    T_conc


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
    self.assertTrue(abs(self.model._rr.R - 0.3e-6) <= 1e-9)
    self.model.run()


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

  def testBasic(self):
    receptor = ReceptorStates(self.model.getResult())
    self.assertEqual(len(receptor._states), 12)
    # Select a state that should have no occupancy
    state = [s for s in receptor._states if
                s.is_ligand and (s.methylation == 2) and
                s.is_phosphorylated][0]
    b = all([x < 0.01 for x in state.getNominalData()])
    self.assertTrue(b)

  def _testSelectStates(self, func, expected):
    receptor = ReceptorStates(self.model.getResult())
    states = receptor.selectStates(func)
    self.assertEqual(len(states), expected)

  def testSelectStates(self):
    func = (lambda l,p,m: l)
    self._testSelectStates(func, 6)
    func = (lambda l,p,m: p)
    self._testSelectStates(func, 6)
    func = (lambda l,p,m: m == 2)
    self._testSelectStates(func, 4)

  # TODO: Losing about 3% of the receptors
  def testSumStates(self):
    # Test that the sum of all states results in the Receptor conc
    func = (lambda l,p,m: True)
    receptor = ReceptorStates(self.model.getResult())
    total = receptor.sumStates(func)
    T = float(T_conc)
    b = all([abs(T-x)/T < 0.04 for x in total])
    self.assertTrue(b)

  def testFrcStates(self):
    # Test that the sum of all states results in the Receptor conc
    func = (lambda l,p,m: True)
    receptor = ReceptorStates(self.model.getResult())
    total = receptor.frcStates(func)
    b = all([abs(1-x) < 0.01 for x in total])
    self.assertTrue(b)
    


if __name__ == '__main__':
  unittest.main()
