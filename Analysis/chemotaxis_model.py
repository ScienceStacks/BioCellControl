"""
Model of E. coli chemotaxis from Spiro, PNAS, 1997.
Assumptions:
1. Ln, LTn have the same kinetic constants for dephosphorylation
2. Start with all T as T2
3. Dephosphorylation occurs at the same rate for bound and unbound receptors
4. The system begins with all Tar in state T2 (doubly methylated)
5. k1a is estimated based on the rates for k_on for the
   bimolecular reaction L + T -> LT, k_on = k5 = 7*10e7.
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
import sys
sys.path.append("../SbStar/Code")
import tellurium as te
from sbstar import SbStar

SIM_START = 0  # Simulation start time
SIM_END = 500  # Simulation end time
SIM_SAMPLES_PER_TIME = 10
TIME = "time"
TEMPLATE = "chemotaxis.tmpl"
# Concentration constants
B_CONC = "1.7e-6"
R_CONC = "0.3e-6"
T2_CONC = "8e-6"
Y_CONC = "20e-6"
Z_CONC = "40e-6"



class ChemotaxisModel(object):

  MODEL_START = "model chemotaxis\n"
  MODEL_END = "end\n"

  def __init__(self, template=TEMPLATE):
    self._rr = None
    self._template = template
    self._template_model = None
    self._antimony_model = None
    self._setupModel()
    self._result = None
    self._factory = None  # State aggregation factory

  def _setupModel(self):
    """
    Processes the template model
    reads: _template
    writes: _template_model, _antimony_model
    """
    lines = []
    with open(self._template) as file:
      lines.append(file.read())
    self._template_model = "\n".join(lines)
    names = ['B', 'R', 'T2', 'Y', 'Z']
    values = [B_CONC, R_CONC, T2_CONC, Y_CONC, Z_CONC]
    pairs = zip(names, values)
    suffixes = ["%s = %s" % (n, v) for n, v in pairs]
    for suffix in suffixes:
      self._template_model += suffix + "\n"
    sbstar = SbStar(self._template_model)
    self._antimony_model = sbstar.expand()

  def _assembleModel(self):
    """
    Assembles the final the model
    """
    return "%s\n%s\n%s" % (ChemotaxisModel.MODEL_START,  \
        self._antimony_model, ChemotaxisModel.MODEL_END)

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
    self._antimony_model = "%s\n%s"  \
        % (self._antimony_model, stg)

  def getModel(self, is_template=True):
    """
    :param bool is_template: return the templated model
    :return str:
    """
    if is_template:
      return self._template_model
    else:
      return self._antimony_model

  def getReceptorStates(self):
    return ReceptorStates(self._result)

  def getYpFraction(self):
    var_Yp = self.getVariable('Yp')
    var_Y = self.getVariable('Y')
    total_Y = var_Yp + var_Y
    return var_Yp/total_Y

  def getBpFraction(self):
    var_Bp = self.getVariable('Bp')
    var_B = self.getVariable('B')
    total_B = var_Bp + var_B
    return var_Bp/total_B

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

  def __init__(self, is_boundL, is_phosphorylated, is_boundR, \
      methylation, simulation_result):
    """
    :param bool is_boundL: bound to a ligand
    :param bool is_phosphorylated: is phosphorylated
    :param bool is_boundR: bound to a CheR
    :param int methylation: methylation level
    :param StructuredArray simulation_result:
    """
    self.is_boundL = is_boundL
    self.is_phosphorylated = is_phosphorylated
    self.is_boundR = is_boundR
    self.methylation = methylation
    self._name = self._makeName()
    self._data = self._extractData(simulation_result)

  def _makeName(self):
    if self.is_boundL:
      name = "LT"
    else:
      name = "T"
    name = "%s%d" % (name, self.methylation)
    if self.is_phosphorylated:
      name = name + "p"
    if self.is_boundR:
      name = name + "R"
    #return "[" + name + "]"
    return name

  def _extractData(self, simulation_result):
    return simulation_result[self._name]

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
  The state has 4 components:
    Is bound to a Ligand (boolean)
    Is phosphorylated (boolean)
    Is bound to CheR (boolean)
    Methylation level (int)
  """

  def __init__(self, result):
    """
    Takes as input a structured array from a simulation result.
    """
    self._result = result
    self._makeStates()

  def _makeStates(self):
    is_boundLs = [False, True]
    is_phosphorylateds = [False, True]
    is_boundRs = [False, True]
    methylations = [2, 3, 4]
    states = []
    for l in is_boundLs:
      for p in is_phosphorylateds:
        for r in is_boundRs:
          for m in methylations:
            states.append(State(l, p, r, m, self._result))
    self._states = states
    total = self.sumStates(lambda l, p, r, m: True)
    [s.setNominalValue(total) for s in states]

  def selectStates(self, func):
    """
    Returns a list of states for which func is true
    :param bool-Function func: arguments are is_boundL, is_phosphorylated, methylation
    :return list-of-State: those states for which func returns True
    """
    return [s for s in self._states
            if func(s.is_boundL, s.is_phosphorylated, s.is_boundR, s.methylation)]

  def sumStates(self, func):
    """
    Returns a the sum of values of selected states
    :param bool-Function func: arguments are is_boundL, is_phosphorylated, methylation
    :return ndarray:
    """
    data = [s.getData() for s in self.selectStates(func)]
    return sum(data)

  def frcStates(self, func):
    """
    Returns the fraction of time in the specified states
    :param bool-Function func: arguments are is_boundL, is_phosphorylated, methylation
    :return ndarray:
    """
    data = [s.getNominalData() for s in self.selectStates(func)]
    return sum(data)

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
    is_valid = False
    if letter ==  "_":
      func = (lambda x: True)
      is_valid = True
    else:
      try:
        is_valid = validation_function(letter)
      except Exception:
        is_valid = False
      if is_valid:
        if letter == "T":
          func = (lambda x: x)
        elif letter == "F":
          func = (lambda x: not x)
        elif isinstance(int(letter), int):
          func = eval("(lambda x: x == %s)" % letter)
        else:
          is_valid = False
    if not is_valid:
      msg = "In name %s, the character %s in position %d is invalid."  \
          % (name, letter, position)
      raise ValueError(msg)
    return func

  def v(self, name):
    """
    Returns the values for the name
    :param list-of-5-char name:
    """
    # Construct the globals algorithmically using eval
    valid_boolean_function = lambda x: x in ["T", "F"]
    if len(name) != 5:
      raise ValueError("Name must be 5 characters.")
    lfunc = self._getFunc(name, 1, valid_boolean_function)
    pfunc = self._getFunc(name, 2, valid_boolean_function)
    rfunc = self._getFunc(name, 3, valid_boolean_function)
    mfunc = self._getFunc(name, 4,  \
        lambda x: isinstance(int(x), int) and int(x) > 1 and int(x) < 5)
    func = (lambda l, p, r, m: lfunc(l)  \
        and pfunc(p) and rfunc(p) and mfunc(m))
    if name[0] == "t":
      result = self._receptor_states.sumStates(func)
    else:
      result = self._receptor_states.frcStates(func)
    return result
