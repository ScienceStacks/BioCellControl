"""
Model of E. coli chemotaxis from Spiro, PNAS, 1997.
Assumptions:
1. Ln, LTn have the same kinetic constants for dephosphorylation
2. Start with all T as T2
3. Dephosphorylation occurs at the same rate for bound and unbound receptors
4. The system begins with all Tar in state T2 (doubly methylated)
5. k1a is estimated based on the rates for k_on for the bimolecular reaction L + T -> LT, k_on = k5 = 7*10e7.
   We use the constant ktuning to adjust the fraction of k5 that's used
6. Methylation of LTn has the same kinetics as methylation of Tn (not explicit in Spiro)

The classes are:
  ChemotaxisModel provides an environment for doing computational experiments on the
    chemotaxis system. The workflow is:
    model = ChemotaxisModel() - creates an instance with the model
    model.setup() - creates the roadrunner simulation object
    model.run() - does a simulation
    getReceptorStates(), getYpFrac(), getBpFrac() - provides access to data
 ReceptorStates provides information about the state of the receptors in terms of
    ligand binding, phosphorylation, and methylation
 State represents a single element of the state of receptors
"""

import tellurium as te

SIM_START = 0  # Simulation start time
SIM_END = 500  # Simulation end time
SIM_SAMPLES_PER_TIME = 10


class ChemotaxisModel(object):

  MODEL_START = "      model chemotaxis\n"
  MODEL_END = "      end\n"
  MODEL = '''
        J0: $x0 -> L; k0*x0
        # REACTIONS from Table 3
        # Methylation
        #J1: T2R -> T3 + R; k1c*T2R
        #J2: T3R -> T4 + R; k2c*T3R
        #J3: LT2R -> LT3 + R; k3c*LT2R
        #J4: LT3R -> LT4 + R; k4c*LT3R
        J1: T2 + R  -> T3 + R; k1c*R
        J2: T3 + R -> T4 + R; k2c*R
        J3: LT2 + R -> LT3 + R; k3c*R
        J4: LT3 + R -> LT4 + R; k4c*R
        #   Details of Tn + R <-> TnR
        #J5: T2R -> T2 + R; k1a*T2R
        #J6: T2 + R -> T2R; k1b*T2*R
        #J7: T3R -> T3 + R; k2a*T3R
        #J8: T3 + R -> T3R; k2b*T3*R
        #J9: T4R -> T4 + R; k3a*T4R
        #J10: T4 + R -> T4R; k3b*T4*R
        #   Details of LTn + R <-> LTnR (not in Spiro)
        #J11: LT2R -> LT2 + R; k1a*LT2R
        #J12: LT2 + R -> LT2R; k1b*LT2*R
        #J13: LT3R -> LT3 + R; k2a*LT3R
        #J14: LT3 + R -> LT3R; k2b*LT3*R
        #J15: LT4R -> LT4 + R; k3a*LT4R
        #J16: LT4 + R -> LT4R; k3b*LT4*R
        # Demethylation reactions
        J17: T3 + Bp -> T2 + Bp; k_1*T3*Bp
        J18: T4 + Bp -> T3 + Bp; k_2*T4*Bp  
        J19: LT3 + Bp -> LT2 + Bp; k_3*LT3*Bp
        J20: LT4 + Bp -> LT3 + Bp; k_4*LT4*Bp
        # Ligand binding gives details to L + T -> LT, LT -> L + T
        J21: L + T2 -> LT2; k5*T2*L
        J22: L + T3 -> LT3; k5*T3*L
        J23: L + T4 -> LT4; k5*T4*L
        J24: LT2 -> T2 + L; k_5*LT2
        J25: LT3 -> T3 + L; k_5*LT3
        J26: LT4 -> T4 + L; k_5*LT4
        # Autophosphorylation reactions
        J27: T2 -> T2p; k8*T2
        J28: T3 -> T3p; k9*T3
        J29: T4 -> T4p; k10*T4
        J30: LT2 -> LT2p; k11*LT2
        J31: LT3 -> LT3p; k12*LT3
        J32: LT4 -> LT4p; k13*LT4
        # B phosphorylations
        T2p + B -> T2 + Bp; kb*T2p*B
        T3p + B -> T3 + Bp; kb*T3p*B
        T4p + B -> T4 + Bp; kb*T4p*B
        LT2p + B -> LT2 + Bp; kb*LT2p*B
        LT3p + B -> LT3 + Bp; kb*LT3p*B
        LT4p + B -> LT4 + Bp; kb*LT4p*B
        # Y phosphorylations
        T2p + Y -> T2 + Yp; ky*T2p*Y
        T3p + Y -> T3 + Yp; ky*T3p*Y
        T4p + Y -> T4 + Yp; ky*T4p*Y
        LT2p + Y -> LT2 + Yp; ky*LT2p*Y
        LT3p + Y -> LT3 + Yp; ky*LT3p*Y
        LT4p + Y -> LT4 + Yp; ky*LT4p*Y
        # B & Y dephosphorylations
   JBp: Bp -> B; k_b*Bp
        Yp + Z -> Y + Z; k_y*Yp*Z    
        # CONSTANTS from Table 3, except k1a (which is noted above).
        k0 = 0
        ktuning = 1
        k1b = ktuning*k5
        k1a = (1.7e-6)/k1b
        k2a = k1a
        k2b = k1b
        k3a = k1a
        k3b = k1b
        k1c = 0.17
        k2c = 0.1*k1c
        k3c = 30*k2c
        k4c = k3c
        k5  = 7e7
        k_5 = 70
        k_1 = 4e5
        k_2 = 3e4
        k_3 = k_1
        k_4 = k_2
        k8 = 15
        k9 = 3*k8
        k10 = 3.2*k8
        k11 = 0
        k12 = 1.1*k8
        k13 = 0.72*k10
        kb = 8e5
        ky = 3e7
        k_b = 0.35
        k_y = 5e5
        # INITIAL VALUES from Table 2
        B = 1.7e-6
        Bp = 0
        LT2 = 0
        LT2R = 0
        LT3 = 0
        LT3R = 0
        LT4 = 0
        T2p = 0
        T3p = 0
        T4p = 0
        LT2p = 0
        LT3p = 0
        LT4p = 0
        R = 0.3e-6
        T2 = 8e-6
        T3 = 0
        T4 = 0
        T2R = 0
        T3R = 0
        x0 = 1
        Y = 20e-6
        Yp = 0
        Z = 40e-6
  '''

  def __init__(self):
    self._rr = None
    self._model = ChemotaxisModel.MODEL
    self._result = None

  def _assembleModel(self):
    """
    Assembles the final the model
    """
    return "%s\n%s\n%s" % (ChemotaxisModel.MODEL_START, 
        self._model, ChemotaxisModel.MODEL_END)

  def initialize(self):
    self._rr = te.loada(self._assembleModel())
    return self._rr

  def run(self, start=SIM_START, end=SIM_END, samples=None):
    """
    Runs the simulation and saves the results.
    :param int start: simulation start time
    :param int end: simulation end time
    :param int samples: number of samples
    :return StructuredArray:
    """
    if samples is None:
      samples = SIM_SAMPLES_PER_TIME*(end-start)
    self._result = self._rr.simulate(start, end, samples)
    return self._result

  def appendToModel(self, stg):
    """
    Adds the string to the end of the model
    """
    self._model = "%s\n%s" % (self._model, stg)

  def getModel(self):
    return self._model

  def getReceptorStates(self):
    return ReceptorStates(self._result)

  def getYpFraction(self):
    total_Y = self._result['Y'] + self._['Yp']
    return self._result['Yp']/total_Y

  def getBpFraction(self):
    total_B = self._result['B'] + self._['Bp']
    return self._result['Bp']/total_B

  def getResult(self):
    return self._result

  def getReactionRateForId(self, id):
    """
    Provides the reaction rate for the id
    :param str id: name of the reaction
    :return float:
    """
    idx = self._rr.getReactionIds().index(id)
    return self._rr.getReactionRates()[idx]

  def getConcentrationForId(self, id):
    """
    Provides the concentration for the chemical species
    :param str id: 
    :return float:
    """
    idx = self._rr.getFloatingSpeciesIds().index(id)
    return self._rr.getFloatingSpeciesConcentrations()[idx]


class State(object):
  """
  Represents a single element of the Receptor state.
  """

  def __init__(self, is_ligand, is_phosphorylated, methylation,
      simulation_result):
    """
    :param bool is_ligand: bound to a ligand
    :param bool is_phosphorylated: is phosphorylated
    :param int methylation: methylation level
    :param StructuredArray simulation_result:
    """
    self.is_ligand = is_ligand
    self.is_phosphorylated = is_phosphorylated
    self.methylation = methylation
    self._name = self._makeName()
    self._data = self._extractData(simulation_result)

  def _makeName(self):
    if self.is_ligand:
      name = "LT"
    else:
      name = "T"
    name = "%s%d" % (name, self.methylation)
    if self.is_phosphorylated:
      name = name + "p"
    return "[" + name + "]"

  def _extractData(self, simulation_result):
    try:
      return simulation_result[self._name]
    except Exception as e:
      import pdb; pdb.set_trace()
      pass

  def getName(self):
    return self._name

  def getData(self):
    return self._data


class ReceptorStates(object):
  """
  Provides information about the state of receptors over time.
  The state has 3 components:
    Is bound to a Ligand (True) or not (False)
    Is phosphorylated (True) or not (False)
    Methylation level (an int)
  """
  
  def __init__(self, result):
    """
    Takes as input a structured array from a simulation result.
    """
    self._result = result
    self._states = self._makeStates()

  def _makeStates(self):
    is_ligands = [False, True]
    is_phosphorylateds = [False, True]
    methylations = [2, 3, 4]
    states = [ [ [State(l, p, m, self._result) for l in is_ligands]
               for p in is_phorphorylateds] for m in methylations]
    return states

  def selectStates(self, func):
    """
    Returns a list of states for which func is true
    :param bool-Function func: arguments are is_ligand, is_phosphorylated, methylation
    :return list-of-State: those states for which func returns True
    """
    return [s for s in self._states 
            if func(s.is_ligand, s.is_phosphorylated, s.methylation)]

  def sumStates(self, func):
    """
    Returns a the sum of values of selected states
    :param bool-Function func: arguments are is_ligand, is_phosphorylated, methylation
    :return ndarray:
    """
    data = [s.getData() for s in self.selectStates(func)]
    return(sum(data))
