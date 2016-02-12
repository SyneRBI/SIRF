import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')
import pGadgetron
import pGadgets

#msec = int(round(time.time() * 1000))
#acq_file = 'acq' + str(msec) + '.h5' 
#print(acq_file)

try:
    acq_proc = pGadgetron.AcquisitionsProcessor()

    print(acq_proc.acq_file)

except pGadgetron.error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
