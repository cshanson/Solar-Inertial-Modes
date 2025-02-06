import os
import errno
import datetime 

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
