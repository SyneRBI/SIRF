import pygadgetron

class error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def _check_status(handle):
    if pygadgetron.executionStatus(handle) != 0:
        msg = pygadgetron.executionError(handle)
        file = pygadgetron.executionErrorFile(handle)
        line = pygadgetron.executionErrorLine(handle)
        errorMsg = \
            repr(msg) + ' exception thrown at line ' + \
            repr(line) + ' of ' + file
        raise error(errorMsg)

class GadgetWithProperties:
    def set_property(self, prop, value):
        handle = pygadgetron.cGT_setGadgetProperty(self.handle, prop, value)
        _check_status(handle)
        pygadgetron.deleteDataHandle(handle)

class GadgetIsmrmrdAcquisitionMessageReader:
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject\
            ('GadgetIsmrmrdAcquisitionMessageReader')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class MRIImageReader:
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('MRIImageReader')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class MRIImageWriter:
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('MRIImageWriter')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class RemoveROOversamplingGadget:
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('RemoveROOversamplingGadget')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class AcquisitionAccumulateTriggerGadget(GadgetWithProperties):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject\
            ('AcquisitionAccumulateTriggerGadget')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class BucketToBufferGadget(GadgetWithProperties):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('BucketToBufferGadget')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class SimpleReconGadget:
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('SimpleReconGadget')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class ImageArraySplitGadget:
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('ImageArraySplitGadget')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class ExtractGadget:
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('ExtractGadget')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class ImageFinishGadget:
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('ImageFinishGadget')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

class SimpleReconGadgetSet(GadgetWithProperties):
    def __init__(self):
        self.handle = None
        self.handle = pygadgetron.cGT_newObject('SimpleReconGadgetSet')
        _check_status(self.handle)
    def __del__(self):
        if self.handle is not None:
            pygadgetron.deleteObject(self.handle)

