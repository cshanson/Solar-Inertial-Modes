from   scipy.integrate   import quad  as quad
from   scipy.integrate   import simps as simps
from   matplotlib.pyplot import plot  as PLOT
from   .          import *
import numpy as NP

#########################################################################
# Utilities

def oddFactorial(k):
    f = k
    while k >= 3:
        k -= 2
        f *= k
    return f

def oddOnEvenFactorial(k):

  f = 1
  for m in range(1,k+1):
    f *= (2*m+1)/(2.0*m)
  return f

#########################################################################
# Basic legendre polynomial computations 
# It is not recommended to use them in a loop in l, l/m

def legendre(l,x,normalized=True,getDerivative=False):

  if NP.amin(x) < -1.0 or NP.amax(x) > 1.0:
    raise ValueError('abscissa <-1 or >1')
  
  if normalized:
    if hasattr(x,'shape'):
      p_l = 1/NP.sqrt(2.0)*NP.ones(x.shape)
    else:
      p_l = 1/NP.sqrt(2.0)*NP.ones((1,))
  else:
    if hasattr(x,'shape'):
      p_l = 1.0*NP.ones(x.shape)
    else:
      p_l = 1.0*NP.ones((1,))
  if l == 0:
    if getDerivative:
      return p_l,NP.zeros(p_l.shape)
    else:
      return p_l

  # Term l=1
  if normalized:
    p_ll = x*NP.sqrt(1.5) 
  else:
    p_ll = p_l*x

  if getDerivative:
    if normalized:
      dp_ll = NP.sqrt(1.5)*NP.ones(p_ll.shape)
    else:
      dp_ll = NP.ones(p_ll.shape)

  if l == 1:
    if getDerivative:
      return p_ll,dp_ll
    else:
      return p_ll

  # Recursion
  for L in range(2,l+1):
    result = iterLegendre(L,x,p_ll,p_l,normalized)
    if getDerivative:
      resD = iterLegendreDerivative(L,x,result,p_ll,normalized)    
    p_l  = p_ll
    p_ll = result

  if getDerivative:
    return result,resD
  else:
    return result

def associatedLegendre(l,m,x,normalized=True):

  ''' Returns the associated Legendre polynomial,
      normalized by sqrt((2l+1)/2)*sqrt((l-m)!/(l+m)!) if required

      Negative ms calls the function with positive m.

      Recursion is performed from Pmm towards Plm (as |m|<=l).
  '''

  if NP.amin(x) < -1.0 or NP.amax(x) > 1.0:
    raise ValueError('abscissa <-1 or >1')
  if abs(m)>l:
    raise ValueError('degree m not compatible with order l : m=%d, l=%d'%(m,l))
  
  # Negative ms
  if (m<0):
    if normalized:
      return (-1)**abs(m)*associatedLegendre(l,abs(m),x,normalized)
    else:
      return (-1)**abs(m)/NP.product(NP.arange(l-abs(m)+1.e0,l+abs(m)+1.e0))*associatedLegendre(l,abs(m),x,normalized)

  # Initial term of the recursion
  if m == 0:
    if normalized:
      if hasattr(x,'shape'):
        p_mm = 1/NP.sqrt(2.0)*NP.ones(x.shape)
      else:
        p_mm = 1/NP.sqrt(2.0)*NP.ones((1,))
    else:
      if hasattr(x,'shape'):
        p_mm = 1.e0*NP.ones(x.shape)
      else:
        p_mm = 1.e0*NP.ones((1,))
  else:
    s = 1
    if m & 1:
      s = -1
    z = NP.sqrt(1-x**2)
    if normalized:
      p_mm = s*z**m*NP.sqrt(oddOnEvenFactorial(m)/2.0)
    else:
      p_mm = s*oddFactorial(2*m-1)*z**m

  if m == l:
    return p_mm

  # Next terms
  p_lm  = p_mm
  if normalized:
    p_llm = p_mm*x*NP.sqrt(2.0*m+3) 
  else:
    p_llm = p_mm*x*(2*m+1)

  if l == m+1:
    return p_llm 

  # Recursion
  for L in range(m+2,l+1):
    result = iterLAssociatedLegendre(L,m,x,p_llm,p_lm,normalized)
    p_lm   = p_llm
    p_llm  = result

  return result

def sphericalHarmonic(ls,ms,theta,phi,normalized=True):

  ''' theta and phi can be given as 1D vectors or meshgrid arrays '''
  th = theta
  f  = phi

  if not hasattr(ls,'__len__'):
    ls = NP.array([ls])
  if not hasattr(ms,'__len__'):
    ms = NP.array([ms])
  if hasattr(th,'__len__'):
    if th.ndim > 1:
      th = th[:,0]
  else:
    th = NP.array([th])
  if hasattr(f,'__len__'):
    if f.ndim > 1:
      f = f[0,:]
  else:
    f = NP.array([f])

  # Plm = associatedLegendre(l,m,NP.cos(th),normalized)
  Plm = associatedLegendre_grid(ls,ms,theta,normalized)
  Eim = NP.exp(1j*ms[:,NP.newaxis]*f[NP.newaxis,:])
  if normalized:
    Eim /= NP.sqrt(2.e0*NP.pi)

  return Plm[...,NP.newaxis]*Eim[NP.newaxis,:,NP.newaxis,:]

def sphericalHarmonic_old(l,m,theta,phi,normalized=True):

  ''' theta and phi can be given as 1D vectors or meshgrid arrays '''
  th = theta
  f  = phi

  if hasattr(th,'__len__'):
    if th.ndim > 1:
      th = th[:,0]
  else:
    th = NP.array([th])
  if hasattr(f,'__len__'):
    if f.ndim > 1:
      f = f[0,:]
  else:
    f = NP.array([f])

  Plm = associatedLegendre(l,m,NP.cos(th),normalized)
  Eim = NP.exp(1j*m*f)
  if normalized:
    Eim /= NP.sqrt(2.e0*NP.pi)

  return Plm[:,NP.newaxis]*Eim[NP.newaxis,:]


#########################################################################
# Iterations method to go from Pl to Pl+1, Plm to Pl+1m

def iterLegendre(L,x,PLm1,PLm2,normalized=True):
  ''' given P_(L-1) and P_(L-2), returns P_(L)'''

  l = float(L)
  if normalized:
    result = (2*l-1)/l*NP.sqrt((2*l+1)/(2*l-1))*x*PLm1 - (l-1)/l*NP.sqrt((2*l+1)/(2*l-3.0))*PLm2
  else:
    result = (x*(2*l-1)*PLm1 - (l-1)*PLm2)/l
  return result

def iterLegendreDerivative(L,x,PL,PLm1,normalized=True):
  ''' given P_(L-1) and P_(L), returns P'_(L)'''

  l = float(L)
  if normalized:
    result        = -l*(x*PL-NP.sqrt((2*l+1.0)/(2*l-1.0))*PLm1)/(1.0-x*x)
    result[x== 1] = l*(l+1)/2. * NP.sqrt((2*l+1.e0)/2.e0)
    result[x==-1] =  (-1)**(L+1)*l*(l+1)/2. * NP.sqrt((2*l+1.e0)/2.e0)
  else:
    result        = -l*(x*PL-PLm1)/(1-x*x)
    result[x== 1] = l*(l+1)/2.
    result[x==-1] = (-1)**(L+1)*l*(l+1)/2.
   
  return result

def iterLAssociatedLegendre(L,M,x,PLm1M,PLm2M,normalized=True):

  l = float(L)
  m = float(M)
  if normalized:
    result  = (2*l-1)*NP.sqrt((2*l+1)/(2*l-1.0)*(l-m)/(l+m*1.0))*x*PLm1M
    result -= (l+m-1)*NP.sqrt((2*l+1)/(2*l-3.0)*(l-m)*(l-m-1.0)/((l+m)*(l+m-1.0)))*PLm2M
    result /= l-m*1.0

  else:
    result = ((2*l-1)*x*PLm1M - (l+m-1)*PLm2M)/(l-m)
  return result

def iterMAssociatedLegendre(L,M,x,PLMm1,PLMm2,normalized=True):

  l = float(L)
  m = float(M)
  result = NP.zeros(PLMm1.shape)
  if normalized:
    result[x!=0] = -2*(m-1)/NP.sqrt(1-x[x!=0]**2)*x[x!=0]*PLMm1[x!=0]\
                   -(l+m-1)*(l-m+2)*PLMm2[x!=0]
  else:
    result[x!=0]  = -2*(m-1)/NP.sqrt(1-x[x!=0]**2)*x[x!=0]*PLMm1[x!=0] - PLMm2[x!=0]
    result       *= 1.e0/NP.sqrt((l+m)*(l-m+1.))
  return result

#########################################################################
# Projections

def projectOnLegendre(data,Ls,normalized=False,axisTheta=-1,pgBar=False,theta = None):
  ''' Compute the projection coefficients of data on Legendre polynomials,
      along axis_theta.
      data can be a list of data arrays, 
      in order to regroup polynomials computation.
  '''

  if not isinstance(data,list):
    dataL = [data]
  else:
    dataL = data

  Nth  = dataL[0].shape[axisTheta]
  if theta is None:
    th   = NP.linspace(0,NP.pi,Nth)
  else:
    th   = theta
  PLm2 = legendre(0,NP.cos(th),normalized)
  PLm1 = legendre(1,NP.cos(th),normalized) 

  l_now = 2
  Nl    = len(Ls)
  
  if dataL[0].ndim == 1:
    data2 = []
    res   = []
    for d in dataL:
      res  .append(NP.zeros(Nl,dtype=dataL[0].dtype))
      data2.append(dataL)
  # Put theta in last position
  else:
    data2 = []
    res   = []
    for d in dataL:
      data2.append(NP.rollaxis(d,axisTheta,d.ndim))
      dims     = list(data2[-1].shape)
      dims[-1] = Nl
      res.append(NP.zeros(tuple(dims),dtype=d.dtype))

  if pgBar:
    PB = progressBar(Nl,'serial')

  for i in range(Nl):

    # Get correct legendre polynomial
    L = Ls[i]
    if L==0:
      Lgdr = PLm2
    elif L==1:
      Lgdr = PLm1
    else:
      while(l_now <= L):
        PL     = iterLegendre(l_now,NP.cos(th),PLm1,PLm2,normalized)
        PLm2   = PLm1
        PLm1   = PL
        l_now += 1
      Lgdr = PL

    # Project input on Pl
    for iData in range(len(dataL)):
      res[iData][...,i] = simps(data2[iData]*Lgdr*NP.sin(th),th)\
                         /simps(Lgdr**2*NP.sin(th),th)

    if pgBar:
      PB.update()
  if pgBar:
    del PB

  # Put theta axis back in place
  if dataL[0].ndim != 1:
    for i in range(len(res)): 
      res[i] = NP.rollaxis(res[i],res[i].ndim-1,axisTheta)

  if len(dataL)==1:
    return res[0]
  else:
    return res

def projectOnAssociatedLegendre(data,ls,ms,theta=None,normalized=True,axisTheta=-1,pgBar=False):
  '''Project data onto its associated Legendre coefficients f_{lm} = \int f(\theta) P_l^m(\theta) sin(\theta) d\theta.'''


  # Theta
  Nth  = data.shape[axisTheta]
  if theta is None:
    th = NP.linspace(0,NP.pi,Nth)
  else:
    th = theta

  # Set geometry and modes
  if isinstance(ls,int):
    Ls = [ls]
  else:
    Ls = ls
  if isinstance(ms,int):
    Ms = [ms]
  else:
    Ms = list(ms)
  x = NP.cos(th)

  # Initiate results arrays and roll the axes of inputs
  # to have L as second last dimension and
  # M as last dimension
  data = NP.array(data)
  data = NP.moveaxis(data,axisTheta,-1)

  if pgBar:
    N  = 0 
    for l in Ls:
      N += sum(abs(NP.array(Ms))<=l)
    PB = progressBar(N,'serial')

  dimData = data.shape[:-1]
  dimData = NP.append(NP.array(dimData),(len(Ls),len(Ms)))
  res = NP.zeros(dimData, dtype=data.dtype)

  # Recursion is done with Plm(Plm-1,Plm-2), from Pmm to Plm, 
  # so the outer loop is in m. We get Pl-m at the same time

  for im in range(len(Ms)):

    m  = Ms[im]
    am = abs(m)

    # Get Pmm and Pm+1m to start recursion (l=m)
    if am == 0:
      if normalized:
        pLM = 1/NP.sqrt(2.0)*NP.ones(x.shape)
      else:
        pLM = 1.e0*NP.ones(x.shape)
    else:
      z = NP.sqrt(1-x**2)
      if normalized:
        pLM = (-1)**am*z**am*NP.sqrt(oddOnEvenFactorial(am)/2.0)
      else:
        pLM = (-1)**am*oddFactorial(2*am-1)*z**am

    if normalized:
      pLp1M = pLM*x*NP.sqrt(2.0*am+3) 
    else:
      pLp1M = pLM*x*(2*am+1)

    # Loop in L
    Lnow = am+2
    for il in range(len(Ls)):

      # Compute only if l>=m
      l = Ls[il]

      if l >= am:

        # Get appropriate associated Legendre polynomial
        if l==am:
          ALgdr = pLM
        elif l==am+1:
          ALgdr = pLp1M
        else:
          while(Lnow <= l):
            pLp2M = iterLAssociatedLegendre(Lnow,am,x,pLp1M,pLM,normalized)
            pLM   = pLp1M
            pLp1M = pLp2M
            Lnow += 1
            ALgdr = pLp2M

        if m<0:
          if normalized:
            ALgdr = ALgdr*(-1)**am
          else:
            ALgdr = ALgdr*(-1)**am/(1.e0*NP.product(NP.arange(l-am+1,l+am+1)))

        # Project input on ALgdr
        res[...,il,im] = simps(data*ALgdr*NP.sin(th),th)
        if not normalized:
          tmpres = simps(ALgdr**2*NP.sin(th),th)
          res[...,il,im] /= tmpres

        del ALgdr
        if pgBar:
          PB.update()

  if pgBar:
    del PB
 
  # Revert axes to their previous place
  #res = NP.moveaxis(res,-1,axisTheta)
  return res


def projectOnSphericalHarmonics(data,ls,ms,theta=None,phi=None,normalized=True,axisTheta=-2,axisPhi=-1,pgBar=False,nbCores=1):
  '''Project data onto its spherical harmonic coefficients f_{lm} = \int f(\theta,\phi) Y_l^m(\theta,\phi) sin(\theta) d\theta d\phi.'''
  isParallelizable = (data.ndim > 2 and NP.prod(data.shape[:-2]) > nbCores)
  data = NP.moveaxis(data,[axisTheta,axisPhi],[-2,-1]) # theta,phi are now the last two dimensions
  if nbCores != 1 and isParallelizable: # parallel version
    dataShape = data.shape
    data = NP.reshape(data, (NP.prod(data.shape[:-2]), data.shape[-2], data.shape[-1])) # data is now N,theta,phi where N contains all the extra dimensions
    data = NP.moveaxis(data,0,-1) # data is now theta,phi,N
    dataLM = reduce(projectOnAssociatedLegendre, (data,ls,ms,theta,normalized,0),data.shape[-1], nbCores, None,progressBar=pgBar) # dataLM is phi,l,m,N
    dataLM = NP.moveaxis(dataLM,-1,0) #  dataLM is N,phi,l,m
    dataLM = NP.reshape(dataLM, (list(dataShape[:-2]) + [dataShape[-1], len(ls), len(ms)])) # dataLM is ...,phi,l,m
  else:
    dataLM = projectOnAssociatedLegendre(data,ls,ms,theta,normalized,axisTheta=-2,pgBar=pgBar) # dataLM is ...phi,l,m
  # Phi
  Nphi  = data.shape[-1]
  if phi is None:
    phi = NP.linspace(0,2.*NP.pi,Nphi)

  dataLM = NP.moveaxis(dataLM,axisPhi-2,-1) 
  dataLM = dataLM * NP.exp(-1.j*ms[:,NP.newaxis]*phi[NP.newaxis,:])
  dataLM = simps(dataLM,phi) # f(l,m)
  return dataLM / NP.sqrt(2.*NP.pi)

#########################################################################
# Reconstructions

def sumLegendre(Al,Ls,x,axisL=0,normalized=True,pgBar=False,sumDerivatives=False):
  ''' Transformation from Legendre polynomial coefficients to actual function in x (cos th)
      Ls: len(Ls) = AL.shape[0] : l-modes on which the coefficients are given
      Several AL arrays can be given in a list (not tuple!)
  '''

  if not isinstance(Al,list):
    AL = [Al]
  else:
    AL = Al

  # Put theta as last dimension
  ALL  = []
  res  = []
  resD = []
  for iCoeffs in range(len(AL)):
    if axisL != -1 or axisL != AL[iCoeffs].ndim-1:
      ALL.append(NP.rollaxis(AL[iCoeffs],axisL,AL[iCoeffs].ndim))
    else:
      ALL.append(AL[iCoeffs])

    dims     = list(ALL[-1].shape)
    dims[-1] = len(x)
    res.append(NP.zeros(tuple(dims),dtype=ALL[-1].dtype)) 
    if sumDerivatives:
      resD.append(NP.zeros(tuple(dims),dtype=ALL[-1].dtype)) 

  PLm2  = legendre(0,x,normalized)
  PLm1  = legendre(1,x,normalized)
  l_now = 2
  if pgBar:
    PB = progressBar(len(Ls),'serial')
  for i in range(len(Ls)):
    # Get correct legendre polynomial
    L = Ls[i]
    if L==0:
      Lgdr = PLm2
      if sumDerivatives:
        dPL = NP.zeros(PLm2.shape)
    elif L==1:
      Lgdr = PLm1
      if sumDerivatives:
        dPL = NP.ones(PLm1.shape)
        if normalized:
          dPL = dPL*NP.sqrt(1.5)
    else:
      while(l_now <= L):
        PL     = iterLegendre(l_now,x,PLm1,PLm2,normalized)
        if sumDerivatives:
          dPL  = iterLegendreDerivative(l_now,x,PL,PLm1,normalized)
        PLm2   = PLm1
        PLm1   = PL
        l_now += 1
      Lgdr = PL
      
    # Add contribution l
    for iCoeffs in range(len(AL)):
      res[iCoeffs] += ALL[iCoeffs][...,i,NP.newaxis]*Lgdr
      if sumDerivatives:
        resD[iCoeffs] += ALL[iCoeffs][...,i,NP.newaxis]*dPL

    if pgBar:
      PB.update()
  if pgBar:
    del PB

  # Put back theta axis
  for iCoeffs in range(len(AL)):
    if axisL != -1 or axisL == AL[iCoeffs].ndim-1:
      res[iCoeffs] = NP.rollaxis(res[iCoeffs],res[iCoeffs].ndim-1,axisL)
      if sumDerivatives:
        resD[iCoeffs] = NP.rollaxis(resD[iCoeffs],res[iCoeffs].ndim-1,axisL)

  if len(AL)==1:
    if sumDerivatives:
      return res[0],resD[0]
    else:
      return res[0]
  else:
    if sumDerivatives:
      return res,resD
    else:
      return res



def reconstructFromAssociatedLegendre(data,ls,ms,theta,normalized=True,axisL=-2,axisM=-1,pgBar=False):
  '''Reconstruct data from its associated Legendre decomposition f(theta) = \sum f_{lm} P_l^m(\theta).'''
  return reconstructFromSphericalHarmonics(data,ls,ms,theta,None,normalized,axisL,axisM,pgBar)


def reconstructFromSphericalHarmonics(data,ls,ms,theta,phi,normalized=True,axisL=-2,axisM=-1,pgBar=False,nbCores=1):
  '''Reconstruct data from its spherical harmonic coefficients f(theta,phi) = \sum f_{lm} Y_l^m(\theta,\phi). If phi is None reconstruct from the associated Legendre.'''

  # Set geometry and modes
  if isinstance(ls,int):
    Ls = [ls]
  else:
    Ls = ls
  if isinstance(ms,int):
    Ms = [ms]
  else:
    Ms = list(ms)
  x = NP.cos(theta)

  # Initiate results arrays and roll the axes of inputs
  # to have L as second last dimension and
  # M as last dimension
  data = NP.array(data)
  data = NP.moveaxis(data,[axisL,axisM],[-2,-1])

  dataShape = data.shape[:-2]
  data = NP.reshape(data, (int(NP.prod(data.shape[:-2])), data.shape[-2], data.shape[-1])) # data is now N,l,m where N contains all the extra dimensions

  if (nbCores > 1 and data.shape[0] > nbCores): # parallel
    data = NP.moveaxis(data,0,-1)  # data is now l,m,N
    res  = reduce(reconstructFromSphericalHarmonicsFunc,(data[NP.newaxis,...],Ls,Ms,theta,phi,normalized,False), data.shape[-1],nbCores,None,progressBar=pgBar)
    res  = NP.array(res)[0,...] # res is theta,phi,N
    res  = NP.moveaxis(res,-1,0) # res is N,theta,phi 
  else: # serial 
    # Call the reconstruction function with axisL=-2 and axisM=-1
    res = reconstructFromSphericalHarmonicsFunc(data,Ls,Ms,theta,phi,normalized,pgBar) 


  # Reshape the result and revert axes to their previous place
  res = NP.reshape(res, (list(dataShape) + [len(theta), len(phi)]))  
  res = NP.moveaxis(res,[-2,-1],[axisL,axisM])

  return res * NP.sqrt(2.*NP.pi)


def reconstructFromSphericalHarmonicsFunc(data,Ls,Ms,theta,phi,normalized=True,pgBar=False):

  if phi is None:
    dimData = NP.append(NP.array(data.shape[0]),len(theta))
  else:
    dimData = NP.append(NP.array(data.shape[0]),(len(theta),len(phi)))
  res = NP.zeros(NP.array(dimData).astype(int), dtype='complex')

  if pgBar:
    N  = 0 
    for l in Ls:
      N += sum(abs(NP.array(Ms))<=l)
    PB = progressBar(N,'serial')

  # Recursion is done with Plm(Plm-1,Plm-2), from Pmm to Plm, 
  # so the outer loop is in m. We get Pl-m at the same time
  x = NP.cos(theta)
  for im in range(len(Ms)):

    m  = Ms[im]
    am = abs(m)

    # Get Pmm and Pm+1m to start recursion (l=m)
    if am == 0:
      if normalized:
        pLM = 1/NP.sqrt(2.0)*NP.ones(x.shape)
      else:
        pLM = 1.e0*NP.ones(x.shape)
    else:
      z = NP.sqrt(1-x**2)
      if normalized:
        pLM = (-1)**am*z**am*NP.sqrt(oddOnEvenFactorial(am)/2.0)
      else:
        pLM = (-1)**am*oddFactorial(2*am-1)*z**am

    if normalized:
      pLp1M = pLM*x*NP.sqrt(2.0*am+3) 
    else:
      pLp1M = pLM*x*(2*am+1)

    # Loop in L
    Lnow = am+2
    for il in range(len(Ls)):

      # Compute only if l>=m
      l = Ls[il]

      if l >= am:

        # Get appropriate associated Legendre polynomial
        if l==am:
          ALgdr = pLM
        elif l==am+1:
          ALgdr = pLp1M
        else:
          while(Lnow <= l):
            pLp2M = iterLAssociatedLegendre(Lnow,am,x,pLp1M,pLM,normalized)
            pLM   = pLp1M
            pLp1M = pLp2M
            Lnow += 1
            ALgdr = pLp2M

        if m<0:
          if normalized:
            ALgdr = ALgdr*(-1)**am
          else:
            ALgdr = ALgdr*(-1)**am/(1.e0*NP.product(NP.arange(l-am+1,l+am+1)))

        # Reconstruct input
        if normalized:
          if phi is None:
            res += data[:,il,im]*ALgdr
          else:
            #print data.shape, ALgdr.shape, phi.shape
            res += data[:,il,im][:,NP.newaxis,NP.newaxis]*(ALgdr[:,NP.newaxis]* NP.exp(1j*m*phi[NP.newaxis,:]))[NP.newaxis,...]
        else:
          tmpres = simps(ALgdr**2*NP.sin(theta),theta)
          if phi is None:
            res += data[...,il,im]*ALgdr*tmpres
          else:
            res += data[...,il,im]*ALgdr[:,NP.newaxis]*tmpres[:,NP.newaxis]* NP.exp(1j*im*phi[NP.newaxis,:])

        del ALgdr
        if pgBar:
          PB.update()

  if pgBar:
    del PB

  if phi is not None:
    res = res / (2.*NP.pi) # ratio between spherical harmonics and Plm
  return res


def sphericalConvolution(data,data0,ls,ms,theta,phi,normalized=True,axisTheta=-2,axisPhi=-1,pgBar=False,nbCores=1):
  '''Compute the spherical convolution between data and data0. data is a function of theta and phi and data0 must depend only on theta (to have only m=0).'''
  # Project 
  datalm = projectOnSphericalHarmonics(data,ls,ms,theta,phi,normalized=normalized,axisTheta=axisTheta,axisPhi=axisPhi,pgBar=pgBar,nbCores=nbCores)
  data0 = data0[...,NP.newaxis]
  data0l = projectOnSphericalHarmonics(data0,ls,NP.array([0]),theta,phi,normalized=normalized,axisTheta=axisTheta,axisPhi=axisPhi,pgBar=pgBar,nbCores=nbCores)

  # Multiply in lm space
  data0l = data0l[...,0]
  convlm = datalm * data0l[...,NP.newaxis] * NP.sqrt(4.*NP.pi/(2*ls[:,NP.newaxis]+1.))

  # Reconstruct
  convtp  = reconstructFromSphericalHarmonics(convlm,ls,ms,theta,phi,normalized=normalized,pgBar=pgBar,nbCores=nbCores) / NP.sqrt(2.*NP.pi)

  return convtp


def legendreArray(n,theta,normalization = False):
  x=NP.cos(theta)
  P=NP.zeros((n+1,len(theta)))
  dP=NP.zeros((n+1,len(theta)))
  if normalization is True:
   P[0,:]=NP.ones((len(x)))/NP.sqrt(2)
   dP[0,:]=NP.zeros((len(x)))
   P[1,:]=x*NP.sqrt(3.0/2)
   dP[1,:]=-NP.sin(theta)*NP.sqrt(3.0/4)

   for l in range(1,n):
     P[l+1,:]=NP.sqrt(2*l+3.)*(NP.sqrt(2*l+1.)*x*P[l,:]-l/NP.sqrt(2*l-1.)*P[l-1,:])/(l+1)
   for l in range(2,n+1):
     dP[l,:]=-NP.sin(theta)*(NP.sqrt(l*(2*l+1.)/(l+1.)))*(-x*P[l,:]/NP.sqrt(2*l+1)+P[l-1,:]/NP.sqrt(2*l-1))/(1-x**2)
     for j in range(0,len(x)):
       if abs(x[j])==1:
         dP[l,j]=0

  elif normalization is False:
    P[0,:]=NP.ones((len(x)))
    dP[0,:]=NP.zeros((len(x)))
    P[1,:]=x
    dP[1,:]=-NP.sin(theta)
    for l in range(1,n):
      P[l+1,:]=((2*l+1.)*x*P[l,:]-l*P[l-1,:])/(l+1.)
    for l in range(2,n+1):
      dP[l,:]=-NP.sin(theta)*(-l*x*P[l,:]+l*P[l-1,:])/(1.-x**2)
      for j in range(0,len(x)):
        if abs(x[j])==1:
          dP[l,j]=0;
          
     
  return P, dP


def Gaunt_Coeffs(l1,l2,l3,m1,m2,m3,Option = None,shift=0):
  Gaunt_exe = pathToMPS() + '/bin/HelioCluster/gaunt_v2.x'
  if Option == 'Plm_Ktheta':
    subP    = subprocess.Popen('%s %s %i %i %i %i' % (Gaunt_exe,Option,l1,l2,l3,shift),shell=True,\
             stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  else:
    subP    = subprocess.Popen('%s %s %i %i %i %i %i %i' % (Gaunt_exe,Option,l1,l2,l3,m1,m2,m3),shell=True,\
               stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  out,err = subP.communicate()
  return NP.array(out.split()).astype(float)

def associatedLegendre_grid(ls,ms,theta=None,normalized=True,pgBar=False):
  '''Project data onto its associated Legendre coefficients f_{lm} = \int f(\theta) P_l^m(\theta) sin(\theta) d\theta.'''


  # Theta
  if theta is None:
    th = NP.linspace(0,NP.pi,1000)
  else:
    th = theta

  # Set geometry and modes
  if isinstance(ls,int):
    Ls = [ls]
  else:
    Ls = ls
  if isinstance(ms,int):
    Ms = [ms]
  else:
    Ms = list(ms)
  x = NP.cos(th)

  if pgBar:
    N  = 0 
    for l in Ls:
      N += sum(abs(NP.array(Ms))<=l)
    PB = progressBar(N,'serial')

  # Recursion is done with Plm(Plm-1,Plm-2), from Pmm to Plm, 
  # so the outer loop is in m. We get Pl-m at the same time

  res = []

  for im in range(len(Ms)):

    m  = Ms[im]
    am = abs(m)

    # Get Pmm and Pm+1m to start recursion (l=m)
    if am == 0:
      if normalized:
        pLM = 1/NP.sqrt(2.0)*NP.ones(x.shape)
      else:
        pLM = 1.e0*NP.ones(x.shape)
    else:
      z = NP.sqrt(1-x**2)
      if normalized:
        pLM = (-1)**am*z**am*NP.sqrt(oddOnEvenFactorial(am)/2.0)
      else:
        pLM = (-1)**am*oddFactorial(2*am-1)*z**am

    if normalized:
      pLp1M = pLM*x*NP.sqrt(2.0*am+3) 
    else:
      pLp1M = pLM*x*(2*am+1)

    # Loop in L
    Lnow = am+2
    resm = []
    for il in range(len(Ls)):

      # Compute only if l>=m
      l = Ls[il]

      if l >= am:

        # Get appropriate associated Legendre polynomial
        if l==am:
          ALgdr = pLM
        elif l==am+1:
          ALgdr = pLp1M
        else:
          while(Lnow <= l):
            pLp2M = iterLAssociatedLegendre(Lnow,am,x,pLp1M,pLM,normalized)
            pLM   = pLp1M
            pLp1M = pLp2M
            Lnow += 1
            ALgdr = pLp2M

        if m<0:
          if normalized:
            ALgdr = ALgdr*(-1)**am
          else:
            ALgdr = ALgdr*(-1)**am/(1.e0*NP.product(NP.arange(l-am+1,l+am+1)))
        resm.append(ALgdr)
        del ALgdr
        if pgBar:
          PB.update()
      else:
        resm.append(NP.zeros(th.shape)*NP.nan)
    res.append(NP.array(resm))

  if pgBar:
    del PB
 
  # Revert axes to their previous place
  #res = NP.moveaxis(res,-1,axisTheta)
  return NP.swapaxes(NP.array(res),0,1)
