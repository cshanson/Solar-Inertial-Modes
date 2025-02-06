import os
import errno
import datetime 
import numpy as np

RSUN = 696.e6

dashLine = "\n------------------------------------------\n"

def mkdir_p(path):
    ''' Equivalent of mkdir -p '''
    try:
      os.makedirs(path)
    except OSError as exc:  # Python >2.5
      if exc.errno == errno.EEXIST and os.path.isdir(path):
        pass
      else:
        raise


def dt_to_dec(dt):
    """Convert a datetime to decimal year."""
    year_start = datetime.datetime(dt.year, 1, 1)
    year_end = year_start.replace(year=dt.year+1)
    return dt.year + ((dt - year_start).total_seconds() /  # seconds so far
        float((year_end - year_start).total_seconds()))  # seconds in year


#==============================================================================
# Sample vectors
def subSampleVector(vector, subSampling):
  if subSampling > 1:
    newSizeR  = int(np.floor(len(vector) / subSampling))
    newVector = np.zeros(newSizeR)
    for i in range(0, newSizeR):
      newVector[i] = (vector[i*subSampling] + vector[i*subSampling+1])/2

    # Why not return vector[::int(subSampling)] ??

  elif subSampling == 1:
    newVector = vector
  else:
    cpt = 0
    sampling  = int(np.round(1. / subSampling))
    newVector = np.zeros(((len(vector) -1) * sampling + 1))
    for i in range(0, len(vector)-1):
      for j in range(0, sampling):
        newVector[cpt] = vector[i] + j * subSampling * (vector[i+1] - vector[i])
        cpt = cpt + 1
    newVector[cpt] = vector[-1]
  return newVector