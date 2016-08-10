import time
import logging

timestr = time.strftime("%Y%m%d-%H%M%S")
logging.basicConfig(filename='ngs_'+timestr+'.log',
    level=logging.DEBUG, format='%(asctime)s %(message)s')

# Am I doing this right?
import filter
import multimer
import sites