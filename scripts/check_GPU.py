import sys

GLOBAL_PATH='/space1/jh7x3/DNSS_release/';

sys.path.insert(0, GLOBAL_PATH+'libs/cudamat-rbm')

try:
  import cudamat
  print "GPU is detected to run cudamat!";
except ImportError:
  print 'GPU failed to be detected!'


