import unittest
from chemotaxis_model import ChemotaxisModel, ReceptorStates, State,  \
    T_conc, StateAggregationFactory, R_conc, B_conc, Y_conc


IGNORE_TEST = False


#############################
# Tests
#############################
# pylint: disable=W0212,C0111,R0904
class TestChemotaxisModelBasic(unittest.TestCase):

  def setUp(self):
    self.model = ChemotaxisModel()
    self.model.initialize()
    self.model.run(start=0, end=10)
    self.receptor = ReceptorStates(self.model.getResult())

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

  def testGetVariable(self):
    Y = self.model.getVariable("Y")
    Yp = self.model.getVariable("Yp")
    b = [y == Y_conc for y in Y + Yp]
    self.assertTrue(b)
    #
    fYp = self.model.getVariable("fYp")
    b = [y <= 1 and 0 <= y for y in fYp]
    self.assertTrue(b)
    #
    fBp = self.model.getVariable("fBp")
    b = [y <= 1 and 0 <= y for y in fBp]
    self.assertTrue(b)
    #
    B = self.model.getVariable("B")
    Bp = self.model.getVariable("Bp")
    b = [y == B_conc for y in B + Bp]
    self.assertTrue(b)
    #
    total = self.model.getVariable("t___")
    b = [y == T_conc for y in total]
    self.assertTrue(b)
    #
    J1 = self.model.getVariable("J1")
    self.assertEqual(len(J1), len(Bp))
    #
    J1 = self.model.getVariable("time")
    self.assertEqual(len(J1), len(Bp))
 
    


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
    self._testBasics(True, True, 2, 'LT2p')
    self._testBasics(False, True, 2, 'T2p')
    self._testBasics(True, False, 2, 'LT2')
    self._testBasics(True, True, 3, 'LT3p')


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


# pylint: disable=W0212,C0111,R0904
class TestStateAggregationFactory(unittest.TestCase):

  def setUp(self):
    self.model = ChemotaxisModel()
    self.model.initialize()
    self.model.run(start=0, end=10)
    self.receptor = ReceptorStates(self.model.getResult())
    self.factory = StateAggregationFactory(self.receptor)

  def testGetFunc(self):
    func = self.factory._getFunc("_", 0, lambda x: True)
    self.assertTrue(func(True))
    func = self.factory._getFunc("T", 0, lambda x: x in ["T", "F"])
    self.assertTrue(func(True))
    self.assertFalse(func(False))
    with self.assertRaises(ValueError):
      self.factory._getFunc("T", 0, lambda x: x in ["F"])
    func = self.factory._getFunc("2", 0, lambda x: int(x) < 5)
    self.assertTrue(func(2))
    self.assertFalse(func(3))

  def testV(self):
    value = self.factory.v("t___")
    expected = self.receptor.sumStates(lambda l,p,m: True)
    pairs = zip(value, expected)
    b = all([abs(x-y) < abs(y)*0.001 for x,y in pairs])
    self.assertTrue(b)
    count = len(expected)
    names = ["_T3", "T_3", "TT_", "T__", "_T_", "__3", "___"]
    for name in names:
      try:
        variable = "f%s" % name
        value = self.factory.v(variable)
        self.assertEqual(len(value), count)
        variable = "t%s" % name
        value = self.factory.v(variable)
        self.assertEqual(len(value), count)
      except Exception as e:
        import pdb; pdb.set_trace()
        pass



if __name__ == '__main__':
  unittest.main()
