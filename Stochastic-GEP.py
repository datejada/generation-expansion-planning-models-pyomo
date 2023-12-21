# %% Two-Stage Stochastic Generation Expansion Planning

#-------------------------------------------------------------------------
#        Universidad Pontificia Comillas
#        Master's Degree in the Electric Power Industry (MEPI)
#        Diego Alejandro Tejada Arango
#        Fecha: 20/11/2021
#-------------------------------------------------------------------------

#%% Import packages
from __future__ import division
from sys import executable
from pyomo.environ import *
from pyomo.opt import SolverFactory

# %% General parameters

# folders names
input_folder  = 'inputs'
output_folder = 'outputs'

# solver definitions
SolverName     = 'cbc'
SolverPath_exe = 'C:\\cbc-win64\\cbc'

# %% We define an abstract model
mGEP = AbstractModel()

# %% Sets and parameters of the abstract model

# sets
mGEP.p  = Set() #time periods (e.g., hours)
mGEP.sc = Set() #uncertainty scenarios
mGEP.g  = Set() #generation technologies

# parameters
mGEP.pScProb = Param(mGEP.sc) #scenario probability [p.u.]
mGEP.pDemand = Param(mGEP.p ) #demand per time period [MW]

mGEP.pVarCost = Param(mGEP.g) #variable   cost of generation units [kEUR/MWh]
mGEP.pInvCost = Param(mGEP.g) #investment cost of generation units [kEUR/MW/year]
mGEP.pUnitCap = Param(mGEP.g) #capacity        of generation units [MW]
mGEP.pIsRenew = Param(mGEP.g) #renewable units boolean indicator   [0,1]

mGEP.pWeight  = Param(within=NonNegativeReals) #weight of representative period [days]
mGEP.pENSCost = Param(within=NonNegativeReals) #energy not supplied cost    [kEUR/MWh]

mGEP.pAviProf = Param(mGEP.sc,mGEP.g,mGEP.p,default=1,mutable=True)   #availability profile [p.u.]

# %% Variables of the abstract model

mGEP.vInvesCost   = Var(                      domain=NonNegativeReals, doc='Total investment Cost                [kEUR]')
mGEP.vOperaCost   = Var(                      domain=NonNegativeReals, doc='Total operating  Cost                [kEUR]')
mGEP.vProduct     = Var(mGEP.sc,mGEP.g,mGEP.p,domain=NonNegativeReals, doc='generation production per scenario   [MW]  ')
mGEP.vInstalUnits = Var(        mGEP.g,       domain=Integers        , doc='number of installed generation units [N]   ')
mGEP.vENS         = Var(mGEP.sc,       mGEP.p,domain=NonNegativeReals, doc='energy not supplied   per scenario   [MW]  ')

# %% Objective function of the abstract model
def eTotalCost(model):
  return   model.vInvesCost + model.vOperaCost
mGEP.eTotalCost = Objective(rule=eTotalCost)

# %% Constraints

#Total investment Cost [kEUR]
def eInvesCost(model):
  return model.vInvesCost == sum(model.pInvCost[g]*model.pUnitCap[g]*model.vInstalUnits[g] for g in model.g)
mGEP.eInvesCost = Constraint(rule=eInvesCost)

#Total operating  Cost [kEUR]
def eOperaCost(model):
  return model.vOperaCost == (model.pWeight*(sum(model.pScProb[sc]*model.pVarCost[g]*model.vProduct[sc,g,p] for sc,g,p in model.sc*model.g*model.p) +
                                             sum(model.pScProb[sc]*model.pENSCost   *model.vENS    [sc,  p] for sc,  p in model.sc*        model.p)))
mGEP.eOperaCost = Constraint(rule=eOperaCost)

#Power balance constraint [MW]
def eBalance(model,sc,p):
  return sum(model.vProduct[sc,g,p] for g in model.g) + model.vENS[sc,p] == model.pDemand[p]
mGEP.eBalance = Constraint(mGEP.sc,mGEP.p,rule=eBalance)

#Max generation constraint
def eMaxProd(model,sc,g,p):
  return model.vProduct[sc,g,p] <= model.pAviProf[sc,g,p]*model.pUnitCap[g]*model.vInstalUnits[g]
mGEP.eMaxProd = Constraint(mGEP.sc,mGEP.g,mGEP.p,rule=eMaxProd)

#Max ENS constraint
def eENSProd(model,sc,p):
  return model.vENS[sc,p] <= model.pDemand[p]
mGEP.eENSProd = Constraint(mGEP.sc,mGEP.p,rule=eENSProd)

# %% We define the optimization solver. You can also use cplex, gurobi, etc
opt = SolverFactory(SolverName,executable=SolverPath_exe)
opt.options['allowableGap'] = 0 

# %% We open a DataPortal to load the data
data = DataPortal() 

# We read all the data from different files
# Scalars
data.load(filename='.\\'+input_folder+'\\scalars.dat')   

# Sets
data.load(filename='.\\'+input_folder+'\\oGEP_Data_Demand.csv'    ,format='set', set='p' )
data.load(filename='.\\'+input_folder+'\\oGEP_Data_Generation.csv',format='set', set='g' )
data.load(filename='.\\'+input_folder+'\\oGEP_Data_Scenario.csv'  ,format='set', set='sc')

# Parameters
data.load(filename='.\\'+input_folder+'\\oGEP_Data_Demand.csv'    ,index=          'p' , param= 'pDemand' )
data.load(filename='.\\'+input_folder+'\\oGEP_Data_Generation.csv',index=      'g'     , param=['pVarCost','pInvCost','pUnitCap','pIsRenew'])
data.load(filename='.\\'+input_folder+'\\oGEP_Data_Scenario.csv'  ,index= 'sc'         , param= 'pScProb' )
data.load(filename='.\\'+input_folder+'\\oGEP_Data_GenAviProf.csv',index=['sc','g','p'], param= 'pAviProf')

# %% We create an instance  
instance = mGEP.create_instance(data)

# We can display all the info of the instance
instance.pprint()

# %% We solve the optimization problem
results = opt.solve(instance,symbolic_solver_labels=True,tee=True) 

# %% Print results

# Check if the problem is optimal
if results.solver.status == SolverStatus.ok and results.solver.termination_condition == TerminationCondition.optimal:
  # objective function value
  print("total cost: "+str(instance.vInvesCost.value+instance.vOperaCost.value))
  # We write some of the results in a csv file
  f = open('.\\'+output_folder+'\\oGEP_Invest_Result.csv', 'w')
  f.write("g,vInstalUnits,pInstalCap"+"\n")
  for g in instance.g.data():
    f.write(str(g)+","+str(instance.vInstalUnits[g].value)+","+str(instance.pUnitCap[g]*instance.vInstalUnits[g].value)+"\n")
  f.close()
else:
  # Print a message indicating that the problem is not optimal
  print("The problem is not optimal.")
  # Print the solver status
  print("Solver Status: "+str(results.solver.status))
