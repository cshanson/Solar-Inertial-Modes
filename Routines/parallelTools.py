from multiprocess import Pool,Lock,Value
from math            import floor,ceil
import sys
import numpy as np
import itertools as IT

# self explanatory
class progressBar(object):

  def __init__(self,N,stype='serial',bonus=None,w=50):

    if stype == 'parallel':
      self.stype    = 'parallel'
      self.progress = Value('i',0)
    else:
      self.stype    = 'serial'
      self.progress = 0.e0
    self.end      = N
    self.lock     = Lock()

    self.w  = w
    str      = '  0% |'+(w*' ')+'|'
    if bonus is not None:
      self.bonus = bonus
      str = str + ' ' + bonus
    else:
      self.bonus = ''
    self.pbl = len(str)
    print( str,end='\r')

  def update(self,w=50):

    if self.stype=='parallel':
      self.lock.acquire()
      self.progress.value += 1
      val = self.progress.value
    else:
      self.progress += 1
      val = self.progress
    percent = int(floor(val/(1.0*self.end)*100))
    prctstr = '%3d%%' % percent
    nequal  = max(percent*self.w//100-1,0) 
 
    print('\r'+prctstr + ' |' + (nequal* '=') + '|' + ((self.w-nequal-1)*' ') + '|' + self.bonus,end='\r' )
    sys.stdout.flush()

    if self.stype=='parallel':
      self.lock.release()
 
  def __del__(self):

    #print (self.pbl*'\r',)
    #print (3*self.pbl*' ',)
    #print (3*self.pbl*'\r',)
    print()
    sys.stdout.flush()




###############################################################################
# Routines to unroll the list of arguments (pool cannot pass several arguments)

def applyPar(args):
  return apppply(*args)

def runningPar(args):
  return running(*args)

def runningSumPar(args):
  return runningSum(*args)

def runningMinPar(args):
  return runningMin(*args)

def runningMaxPar(args):
  return runningMax(*args)

def runningAvePar(args):
  return runningAve(*args)

###############################################################################
# Process or reduction routines

def apppply(funcname,args):

  if args[0][0]:
    PG = progressBar(len(args),'serial',bonus=' (parallel)')

  for i in range(len(args)):
    funcname(*(args[i][1:]))
    if args[i][0]:
      PG.update()          
 
  if args[i][0]:
    del PG

def running(funcname,args):
  res = []
  if args[0][0]:
    PB = progressBar(len(args),'serial',bonus=' (parallel)')
  for i in range(len(args)):
    res.append(funcname(*(args[i][1:])))
    if args[i][0]:
      PB.update()          
  if args[i][0]:
    del PB
  return res

def runningSum(funcname,args):
  ''' compute the sum from start to end, applying func_name to whatever args'''
  res = 0
  if args[0][0]:
    PG = progressBar(len(args),'serial',bonus=' (parallel)')
  for i in range(len(args)):
    res += funcname(*(args[i][1:]))
    if args[i][0]:
      PG.update()          
  if args[i][0]:
    del PG
  return res

def runningAve(funcname,args):
  ''' compute the sum from start to end, applying func_name to whatever args'''
  res = 0
  if args[0][0]:
    PG = progressBar(len(args),'serial',bonus=' (parallel)')
  for i in range(len(args)):
    res += args[i][1]*funcname(*(args[i][2:]))
    if args[i][0]:
      PG.update()          
  if args[i][0]:
    del PG
  return res

def runningMin(funcname,args):
  res = funcname(*(args[0][1:]))
  if args[0][0]:
    PG = progressBar(len(args),'serial',bonus=' (parallel)')
    PG.update()
  for i in range(1,len(args)):
    res = np.minimum(res,funcname(*(args[i][1:])))
    if args[i][0]:
      PG.update()          
  if args[i][0]:
    del PG
  return res

def runningMax(funcname,args):
  res = funcname(*(args[0][1:]))
  if args[0][0]:
    PG = progressBar(len(args),'serial',bonus=' (parallel)')
    PG.update()
  for i in range(1,len(args)):
    res = np.maximum(res,funcname(*(args[i][1:])))
    if args[i][0]:
      PG.update()          
  if args[i][0]:
    del PG
  return res

########################################################################
# Reduction preparation routines

def buildArgTable(args,N,nbproc,progressBar,weights=None):

  starts = partition(N,nbproc)
  # prepare list of arguments
  arg_table = []
  if weights is not None:
    w = np.genfromtxt(weights)

  for i in range(nbproc):
    arg_table_proc = []    
    for j in range(starts[i],starts[i+1]):

      list_args = []
      # Precise in arguments if progress bar is needed
      if progressBar and i==0:
        list_args.append(True)
      else:
        list_args.append(False)

      # Precise the weights for a weighted average
      if weights is not None:
        list_args.append(w[j])
      for arg in list(args):
        # Detect if the argument is a numpy array with last dimension N
        try:
          if (arg.shape[-1]==N):
            if arg.ndim == 1:
              list_args.append(arg[j])
            else:
              list_args.append(arg[...,j])
          else:
            list_args.append(arg)
        except: 
          list_args.append(arg)
      arg_table_proc.append(tuple(list_args))

    arg_table.append(arg_table_proc)
  return arg_table

def partition(N,Nprocs):

  division = N / float(Nprocs)
  starts   = []
  for i in range(Nprocs):
    starts.append(int(round(division * i)))
  starts.append(min(int(round(division*(i+1))),N))
  return starts

########################################################################
# Top reduction routines

def reduce(funcname,args,N,nbproc,type=None,progressBar=False,SparseStack=None):
  ''' Applies the type of reduction on a pool map of nbprocs
      N : size of array on which we perform computation using the function func_name
      type can be 'SUM' 'MAX' 'MIN' 'AVE'
  '''

  if (type is not None and type[:3] == 'AVE'):
    weights = type[3:]
  else:
    weights = None
  argTable = buildArgTable(args,N,nbproc,progressBar,weights)



  pool = Pool(nbproc)
  print(pool)

  if (type==None):
    resPerProc = pool.map(runningPar,zip(IT.repeat(funcname),argTable))
  elif (type=='SUM'):
    resPerProc = pool.map(runningSumPar,zip(IT.repeat(funcname),argTable))
  elif (type=='MIN'):
    resPerProc = pool.map(runningMinPar,zip(IT.repeat(funcname),argTable))
  elif (type=='MAX'):
    resPerProc = pool.map(runningMaxPar,zip(IT.repeat(funcname),argTable))
  elif (type[:3]=='AVE'):
    resPerProc = pool.map(runningAvePar,zip(IT.repeat(funcname),argTable))

  # Final reduction
  if (type==None):
    if SparseStack is None:
      res2 = []
      for i in range(nbproc):
        for j in range(len(resPerProc[i])):
          res2.append(np.array(resPerProc[i][j]))
      D   = np.array(res2)
      res = np.rollaxis(D,0,D.ndim)
    else:
      res = []
      for i in range(nbproc):
        for j in range(len(resPerProc[i])):
          res.append(resPerProc[i][j])
      if SparseStack.upper() == 'VSTACK':
        return vstack(res)
      elif SparseStack.upper() == 'HSTACK':
        return hstack(res)        

  elif (type=='SUM' or type[:3] == 'AVE'):
    res = np.sum(resPerProc,axis=0)
  elif (type=='MIN'):
    res = np.min(resPerProc,axis=0)
  elif (type=='MAX'):
    res = np.max(resPerProc,axis=0)
  elif (type=='NO'):
    # nothing to return
    print ('ending parallel computation')

  pool.close()
  return np.asarray(res)

def apply(funcname,args,N,nbproc,progressBar=False):
  '''same as reduce, without return '''

  argTable = buildArgTable(args,N,nbproc,progressBar)

  pool = Pool(nbproc)
  pool.apply_async(applyPar,zip(IT.repeat(funcname),argTable))
  pool.close()
