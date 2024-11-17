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
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import highspy

# %% General parameters

# folders names
input_folder  = 'inputs'
output_folder = 'outputs'

# solver definitions
SolverName = 'appsi_highs' #gurobi

# %% We define an abstract model
mGEP = pyo.AbstractModel()

# %% Sets and parameters of the abstract model

# sets
mGEP.p  = pyo.Set() #time periods (e.g., hours)
mGEP.sc = pyo.Set() #uncertainty scenarios
mGEP.g  = pyo.Set() #generation technologies

# parameters
mGEP.pScProb = pyo.Param(mGEP.sc) #scenario probability [p.u.]
mGEP.pDemand = pyo.Param(mGEP.p ) #demand per time period [MW]

mGEP.pVarCost = pyo.Param(mGEP.g) #variable   cost of generation units [kEUR/MWh]
mGEP.pInvCost = pyo.Param(mGEP.g) #investment cost of generation units [kEUR/MW/year]
mGEP.pUnitCap = pyo.Param(mGEP.g) #capacity        of generation units [MW]
mGEP.pIsRenew = pyo.Param(mGEP.g) #renewable units boolean indicator   [0,1]

mGEP.pWeight  = pyo.Param(within=pyo.NonNegativeReals) #weight of representative period [days]
mGEP.pENSCost = pyo.Param(within=pyo.NonNegativeReals) #energy not supplied cost    [kEUR/MWh]

mGEP.pAviProf = pyo.Param(mGEP.sc,mGEP.g,mGEP.p,default=1,mutable=True)   #availability profile [p.u.]

# %% Variables of the abstract model

mGEP.vInvesCost   = pyo.Var(                      domain=pyo.NonNegativeReals, doc='Total investment Cost                [kEUR]')
mGEP.vOperaCost   = pyo.Var(                      domain=pyo.NonNegativeReals, doc='Total operating  Cost                [kEUR]')
mGEP.vProduct     = pyo.Var(mGEP.sc,mGEP.g,mGEP.p,domain=pyo.NonNegativeReals, doc='generation production per scenario   [MW]  ')
mGEP.vInstalUnits = pyo.Var(        mGEP.g,       domain=pyo.Integers        , doc='number of installed generation units [N]   ')
mGEP.vENS         = pyo.Var(mGEP.sc,       mGEP.p,domain=pyo.NonNegativeReals, doc='energy not supplied   per scenario   [MW]  ')

# %% Objective function of the abstract model
def eTotalCost(model):
  return   model.vInvesCost + model.vOperaCost
mGEP.eTotalCost = pyo.Objective(rule=eTotalCost)

# %% Constraints

#Total investment Cost [kEUR]
def eInvesCost(model):
  return model.vInvesCost == sum(model.pInvCost[g]*model.pUnitCap[g]*model.vInstalUnits[g] for g in model.g)
mGEP.eInvesCost = pyo.Constraint(rule=eInvesCost)

#Total operating  Cost [kEUR]
def eOperaCost(model):
  return model.vOperaCost == (model.pWeight*(sum(model.pScProb[sc]*model.pVarCost[g]*model.vProduct[sc,g,p] for sc,g,p in model.sc*model.g*model.p) +
                                             sum(model.pScProb[sc]*model.pENSCost   *model.vENS    [sc,  p] for sc,  p in model.sc*        model.p)))
mGEP.eOperaCost = pyo.Constraint(rule=eOperaCost)

#Power balance constraint [MW]
def eBalance(model,sc,p):
  return sum(model.vProduct[sc,g,p] for g in model.g) + model.vENS[sc,p] == model.pDemand[p]
mGEP.eBalance = pyo.Constraint(mGEP.sc,mGEP.p,rule=eBalance)

#Max generation constraint
def eMaxProd(model,sc,g,p):
  return model.vProduct[sc,g,p] <= model.pAviProf[sc,g,p]*model.pUnitCap[g]*model.vInstalUnits[g]
mGEP.eMaxProd = pyo.Constraint(mGEP.sc,mGEP.g,mGEP.p,rule=eMaxProd)

#Max ENS constraint
def eENSProd(model,sc,p):
  return model.vENS[sc,p] <= model.pDemand[p]
mGEP.eENSProd = pyo.Constraint(mGEP.sc,mGEP.p,rule=eENSProd)

# %% We define the optimization solver. You can also use cplex, gurobi, etc

opt = SolverFactory(SolverName)

# We define the options of the solver (this depends on the solver you are using)
##opt.options['mip_rel_gap'] = 0 # HiGHS option for relative gap

# %% We open a DataPortal to load the data
data = pyo.DataPortal() 

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

# write the optimization problem
instance.write('mGEP.lp', io_options={'symbolic_solver_labels': True})

# %% We solve the optimization problem
results = opt.solve(instance,symbolic_solver_labels=True,tee=True) 

# %% Print results

# Print the number of variables and constraints
print("Number of variables: "+str(instance.nvariables()))
print("Number of constraints: "+str(instance.nconstraints()))

# Check if the problem is optimal
if results.solver.status == pyo.SolverStatus.ok and results.solver.termination_condition == pyo.TerminationCondition.optimal:
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

# %%
