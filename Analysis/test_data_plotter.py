import unittest
from data_plotter import DataPlotter
from chemotaxis_model import ChemotaxisModel


IGNORE_TEST = False


#############################
# Tests
#############################
# pylint: disable=W0212,C0111,R0904
class TestDataPlotter(unittest.TestCase):

  def setUp(self):
    self.model = ChemotaxisModel()
    self.model.initialize()
    self.model.run(start=0, end=10)
    self.plotter = DataPlotter(self.model)

  def testLines(self):
    self.plotter.lines(["fYp", "Bp", "fTT__"])


if __name__ == '__main__':
  unittest.main()
