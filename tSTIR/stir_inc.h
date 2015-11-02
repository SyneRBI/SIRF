#ifndef STIR_DATA_TYPES
#define STIR_DATA_TYPES

#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <numeric>
#include <string>

#include <boost\algorithm\string.hpp>

#include "stir/DiscretisedDensity.h"
#include "stir/ArrayFunction.h"
#include "stir/CartesianCoordinate3D.h"
#include "stir/CPUTimer.h"
#include "stir/DataProcessor.h"
#include "stir/HighResWallClockTimer.h"
#include "stir/IndexRange3D.h"
#include "stir/is_null_ptr.h"
#include "stir/KeyParser.h"
#include "stir/recon_array_functions.h"
#include "stir/Succeeded.h"
#include "stir/utilities.h"
#include "stir/VoxelsOnCartesianGrid.h"
#include "stir/IO/OutputFileFormat.h"
#include "stir/IO/read_from_file.h"
#include "stir/OSMAPOSL/OSMAPOSLReconstruction.h"
#include "stir/OSSPS/OSSPSReconstruction.h"
#include "stir/recon_buildblock/distributable_main.h"
#include "stir/recon_buildblock/PoissonLogLikelihoodWithLinearModelForMeanAndProjData.h"
#include "stir/recon_buildblock/ProjectorByBinPairUsingProjMatrixByBin.h"
#include "stir/recon_buildblock/ProjMatrixByBinUsingRayTracing.h"
#include "stir/recon_buildblock/QuadraticPrior.h"
#include "stir/Shape/Shape3D.h"
#include "stir/TruncateToCylindricalFOVImageProcessor.h"

USING_NAMESPACE_STIR

typedef DiscretisedDensity<3, float> Image3DF;
typedef boost::shared_ptr<Image3DF> ptrImage3DF;
typedef boost::shared_ptr<Image3DF> sptrImage3DF;

#endif