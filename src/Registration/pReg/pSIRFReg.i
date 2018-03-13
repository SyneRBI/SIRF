/* File: pSIRFReg.i */
%module pSIRFReg
%include <std_string.i>
%include <exception.i>

%{
#define SWIG_FILE_WITH_INIT
#include <cmath>
#include "SIRFRegNiftyAladin.h"
#include "SIRFRegNiftyF3d.h"
#include "SIRFRegNiftyResample.h"
#include "SIRFRegImageWeightedMean.h"
#include "SIRFRegActivityCorrect.h"
#include "SIRFRegMisc.h"
%}

%exception {
    try {
        $action
    } catch(const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch(...) {
        SWIG_exception(SWIG_RuntimeError,"Unknown exception");
    }
}

// Misc
%include "../cReg/SIRFRegMisc.h"

// Base class
%include "../cReg/SIRFReg.h"

// Aladin - only offer the user single precision for now
%include "../cReg/SIRFRegNiftyAladin.h"
%template(NiftyAladin) SIRFRegNiftyAladin<float>;
// %template(NiftyAladin_double) SIRFRegNiftyAladin<double>;

// F3d - only offer the user single precision for now
%include "../cReg/SIRFRegNiftyF3d.h"
%template(NiftyF3d) SIRFRegNiftyF3d<float>;
// %template(NiftyF3d_double) SIRFRegNiftyF3d<double>;

// Resample
%rename(NiftyResample) SIRFRegNiftyResample;
%include "../cReg/SIRFRegNiftyResample.h"

// Weighted mean of images
%rename(ImageWeightedMean) SIRFRegImageWeightedMean;
%include "../cReg/SIRFRegImageWeightedMean.h"

// Activity correction
%rename(ActivityCorrect) SIRFRegActivityCorrect;
%include "../cReg/SIRFRegActivityCorrect.h"
