/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

#ifndef STIR_DATA_TYPES
#define STIR_DATA_TYPES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <numeric>
#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

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
#include "stir/ProjDataInfoCylindrical.h"
#include "stir/ProjDataInMemory.h"
#include "stir/ProjDataInterfile.h"
#include "stir/recon_buildblock/BinNormalisationFromProjData.h"
#include "stir/recon_buildblock/distributable_main.h"
#include "stir/recon_buildblock/PoissonLogLikelihoodWithLinearModelForMeanAndProjData.h"
#include "stir/recon_buildblock/ProjectorByBinPairUsingProjMatrixByBin.h"
#include "stir/recon_buildblock/ProjMatrixByBinUsingRayTracing.h"
#include "stir/recon_buildblock/QuadraticPrior.h"
#include "stir/Shape/EllipsoidalCylinder.h"
#include "stir/Shape/Shape3D.h"
#include "stir/TruncateToCylindricalFOVImageProcessor.h"

#include "stir/StirException.h"
#include "stir/TextWriter.h"

#define GRAB 1

USING_NAMESPACE_STIR

typedef DiscretisedDensity<3, float> Image3DF;
typedef boost::shared_ptr<Image3DF> sptrImage3DF;
typedef boost::shared_ptr<ProjData> sptrProjData;
typedef boost::shared_ptr<ProjectorByBinPair> sptrProjPair;
typedef CartesianCoordinate3D<float> Coord3DF;
typedef VoxelsOnCartesianGrid<float> Voxels3DF;
typedef boost::shared_ptr<Voxels3DF> sptrVoxels3DF;
typedef boost::shared_ptr<Shape3D> sptrShape3D;
typedef Reconstruction<Image3DF> Reconstruction3DF;
typedef IterativeReconstruction<Image3DF> IterativeReconstruction3DF;
typedef GeneralisedObjectiveFunction<Image3DF> ObjectiveFunction3DF;
typedef PoissonLogLikelihoodWithLinearModelForMean<Image3DF> 
PoissonLogLhLinModMean3DF;
typedef PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF>
PoissonLogLhLinModMeanProjData3DF;
typedef ProjectorByBinPairUsingProjMatrixByBin ProjectorPairUsingMatrix;
typedef ProjMatrixByBinUsingRayTracing RayTracingMatrix;
typedef GeneralisedPrior<Image3DF> Prior3DF;
typedef QuadraticPrior<float> QuadPrior3DF;
typedef DataProcessor<Image3DF> DataProcessor3DF;
typedef TruncateToCylindricalFOVImageProcessor<float> CylindricFilter3DF;

#endif