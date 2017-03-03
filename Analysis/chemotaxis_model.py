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
TIME = "time"
# Concentration constants
B_conc = "1.7e-6"
R_conc = "0.75e-6"  # Differs from Spiro's 0.3e-6
T_conc = "8e-6"
Y_conc = "20e-6"
Z_conc = "40e-6"



class ChemotaxisModel(object):

  MODEL_START = "      model chemotaxis\n"
  MODEL_END = "      end\n"
  MODEL = '''
        J0: $x0 -> L; k0*x0
        # REACTIONS from Spiro, Table 3
        # Methylation
        J1: T2R -> T3 + R; k1c*T2R
        J2: T3R -> T4 + R; k2c*T3R
        J3: LT2R -> LT3 + R; k3c*LT2R
        J4: LT3R -> LT4 + R; k4c*LT3R
        #   Details of Tn + R <-> TnR
        J5: T2R -> T2 + R; k1a*T2R
        J6: T2 + R -> T2R; k1b*T2*R
        J7: T3R -> T3 + R; k2a*T3R
        J8: T3 + R -> T3R; k2b*T3*R
        #J9: T4R -> T4 + R; k3a*T4R
        #J10: T4 + R -> T4R; k3b*T4*R
        #   Details of LTn + R <-> LTnR (not in Spiro)
        J11: LT2R -> LT2 + R; k1a*LT2R
        J12: LT2 + R -> LT2R; k1b*LT2*R
        J13: LT3R -> LT3 + R; k2a*LT3R
        J14: LT3 + R -> LT3R; k2b*LT3*R
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
        J33:T2p + B -> T2 + Bp; kb*T2p*B
        J34:T3p + B -> T3 + Bp; kb*T3p*B
        J35:T4p + B -> T4 + Bp; kb*T4p*B
        J36:LT2p + B -> LT2 + Bp; kb*LT2p*B
        J37:LT3p + B -> LT3 + Bp; kb*LT3p*B
        J38:LT4p + B -> LT4 + Bp; kb*LT4p*B
        # Y phosphorylations
        J39:T2p + Y -> T2 + Yp; ky*T2p*Y
        J40:T3p + Y -> T3 + Yp; ky*T3p*Y
        J41:T4p + Y -> T4 + Yp; ky*T4p*Y
        J42:LT2p + Y -> LT2 + Yp; ky*LT2p*Y
        J43:LT3p + Y -> LT3 + Yp; ky*LT3p*Y
        J44:LT4p + Y -> LT4 + Yp; ky*LT4p*Y
        # B & Y dephosphorylations
        J45: Bp -> B; k_b*Bp
        J46: Yp + Z -> Y + Z; k_y*Yp*Z    
        # CONSTANTS from Spiro, Table 3, except k1a (which is noted above).
        k0 = 0
        ktuning = 10
        k1b = ktuning*k5
        k1a = (1.7e-6)/k1b
        k2a = k1a
        k2b = k1b
        k3a = k1a
        k3b = k1b
        k1c = 0.17
        k2c = 0.1*k1c
        k3c = 30*k1c
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
        # INITIAL VALUES from Spiro, Table 2
        B = %s
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
        R = %s
        T2 = %s
        T3 = 0
        T4 = 0
        T2R = 0
        T3R = 0
        x0 = 1
        Y = %s
        Yp = 0
        Z = %s
  ''' % (B_conc, R_conc, T_conc, Y_conc, Z_conc)

  def __init__(self):
    self._rr = None
    self._model = ChemotaxisModel.MODEL
    self._result = None
    self._factory = None  # State aggregation factory

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
    selection = self._rr.getFloatingSpeciesIds()
    selection.extend(self._rr.getReactionIds())
    selection.append(TIME)
    self._rr.timeCourseSelections = selection
    self._result = self._rr.simulate(start, end, samples)
    states = self.getReceptorStates()
    self._factory = StateAggregationFactory(states)
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
    Yp = self.getVariable('Yp')
    Y = self.getVariable('Y')
    total_Y = Yp + Y
    return Yp/total_Y

  def getBpFraction(self):
    Bp = self.getVariable('Bp')
    B = self.getVariable('B')
    total_B = Bp + B
    return Bp/total_B

  def getResult(self):
    return self._result

  def _makeVariableName(self, name):
    return "[%s]" % name

  def getVariable(self, name):
    """
    Provides variables in the model.
    :param str name:
    :return np.array:
    """
    result = None
    # Try the special names
    if name == "time":
      result = self._result["time"]
    elif name == "fYp":
      result = self.getYpFraction()
    elif name == "fBp":
      result = self.getBpFraction()
    else:
      # See if it's a state name
      try:
        result = self._factory.v(name)
      except:
        pass
      if result is None:
        # Assume it's a model name
        try:
          if result is None:
            result = self._result[name]
        except:
          pass
    if result is None:
      import pdb; pdb.set_trace()
      raise ValueError("Variable %s not found." % name)
    return result


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
    #return "[" + name + "]"
    return name

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

  def getNominalData(self):
    return self._nominal_data

  def setNominalValue(self, adjustment):
    self._nominal_data = self._data/adjustment


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
    self._makeStates()

  def _makeStates(self):
    is_ligands = [False, True]
    is_phosphorylateds = [False, True]
    methylations = [2, 3, 4]
    states = []
    for l in is_ligands:
      for p in is_phosphorylateds:
        for m in methylations:
          states.append(State(l, p, m, self._result))
    self._states = states
    total = self.sumStates(lambda l,m,p: True)
    [s.setNominalValue(total) for s in states]

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

  def frcStates(self, func):
    """
    Returns the fraction of time in the specified states
    :param bool-Function func: arguments are is_ligand, is_phosphorylated, methylation
    :return ndarray:
    """
    data = [s.getNominalData() for s in self.selectStates(func)]
    return(sum(data))

class StateAggregationFactory(object):
  """
  State aggregations are summations of values of a subset of states.
  An aggregation is
  either a total concentration, denoted by the prefix "t", or a fraction of the total
  concentration of receptors, denoted by "f". States have a three tuple: l, p, m. "l" is
  true (T) if the receptor is bound to a ligand; otherwise it is false. "p" is true if the
  receptor is phosphorylated; otherwise it is false. "m" is the level of methylation. If
  there is no specification for the state component, an underscore is used.
  For example, the fraction of receptors that are bound to ligands and are phosphorylated is
  denoted by "fTT_", and the total concentration of receptors that are bound to ligands,
  and have a methylation of 3 is "tT_3". The total concentration of receptors is all states
  is "t___". Note that "f___" should be 1.
  """

  def __init__(self, receptor_states):
    self._receptor_states = receptor_states

  def _getFunc(self, name, position, validation_function):
    """
    Creates a function to select for a part of a name
    :param str name: name being interpreted
    :param int position: position of the letter in the name
    :param Function validation_function:
    :return Function: function used to select one part of state
    :raises ValueError: invalid character
    """
    letter = name[position]
    if letter ==  "_":
      func = (lambda x: True)
      is_valid = True
    else:
      try:
        is_valid = validation_function(letter)
      except Exception:
        is_valid = False
      if is_valid and letter == "T":
        func = (lambda x: x)
      elif is_valid and letter == "F":
        func = (lambda x: not x)
      elif is_valid and isinstance(int(letter), int):
        func = eval("(lambda x: x == %s)" % letter)
      else:
        is_valid = False
    if not is_valid:
      raise ValueError("In name %s, the character %s in position %d is invalid."
          % (name, letter, position))
    return func
  
  def v(self, name):
    """
    Returns the values for the name
    :param list-of-4-char name:
    """
    # Construct the globals algorithmically using eval
    if len(name) != 4:
      raise ValueError("Name must be 4 characters.")
    lfunc = self._getFunc(name, 1,
        lambda x: x in ["T", "F"])
    pfunc = self._getFunc(name, 2,
        lambda x: x in ["T", "F"])
    mfunc = self._getFunc(name, 3, 
        lambda x: isinstance(int(x),int) and int(x) > 1 and int(x) < 5)
    func = (lambda l,p,m: lfunc(l) and pfunc(p) and mfunc(m))
    if name[0] == "t":
      result = self._receptor_states.sumStates(func)
    else:
      result = self._receptor_states.frcStates(func)
    return result
