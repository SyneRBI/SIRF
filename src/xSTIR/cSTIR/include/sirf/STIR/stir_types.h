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

#include <boost/algorithm/string.hpp>

#include "stir/DataProcessor.h"
#include "stir/DiscretisedDensity.h"
#include "stir/CartesianCoordinate3D.h"
#include "stir/IndexRange3D.h"
#include "stir/is_null_ptr.h"
#include "stir/ML_norm.h"
#include "stir/recon_array_functions.h"
#include "stir/SeparableGaussianImageFilter.h"
#include "stir/Succeeded.h"
#include "stir/utilities.h"
#include "stir/VoxelsOnCartesianGrid.h"
#include "stir/analytic/FBP2D/FBP2DReconstruction.h"
#include "stir/IO/OutputFileFormat.h"
#include "stir/IO/read_from_file.h"
#include "stir/listmode/CListRecord.h"
#include "stir/listmode/CListEventCylindricalScannerWithDiscreteDetectors.h"
#include "stir/listmode/LmToProjData.h"
#include "stir/OSMAPOSL/OSMAPOSLReconstruction.h"
#ifdef USE_HKEM
#include "stir/KOSMAPOSL/KOSMAPOSLReconstruction.h"
#endif
#include "stir/OSSPS/OSSPSReconstruction.h"
#include "stir/ProjDataInfoCylindrical.h"
#include "stir/ProjDataInMemory.h"
#include "stir/ProjDataInterfile.h"
#include "stir/recon_buildblock/BinNormalisationFromAttenuationImage.h"
#include "stir/recon_buildblock/BinNormalisationFromECAT8.h"
#include "stir/recon_buildblock/BinNormalisationFromProjData.h"
#include "stir/recon_buildblock/ChainedBinNormalisation.h"
#include "stir/recon_buildblock/PLSPrior.h"
#include "stir/recon_buildblock/PoissonLogLikelihoodWithLinearModelForMeanAndProjData.h"
#include "stir/recon_buildblock/ProjectorByBinPairUsingProjMatrixByBin.h"
#include "stir/recon_buildblock/ProjMatrixByBinUsingRayTracing.h"
#include "stir/recon_buildblock/QuadraticPrior.h"
#include "stir/Shape/EllipsoidalCylinder.h"
#include "stir/Shape/Shape3D.h"
#include "stir/shared_ptr.h"
#include "stir/SSRB.h"
#include "stir/TruncateToCylindricalFOVImageProcessor.h"

#include "stir/StirException.h"
#include "stir/TextWriter.h"

#define GRAB 1

namespace sirf {

	typedef stir::DiscretisedDensity<3, float> Image3DF;
	typedef stir::shared_ptr<Image3DF> sptrImage3DF;
	typedef stir::shared_ptr<stir::ProjData> sptrProjData;
	typedef stir::CartesianCoordinate3D<float> Coord3DF;
	typedef stir::CartesianCoordinate3D<int> Coord3DI;
	typedef stir::VoxelsOnCartesianGrid<float> Voxels3DF;
	typedef stir::shared_ptr<Voxels3DF> sptrVoxels3DF;
	typedef stir::shared_ptr<stir::Shape3D> sptrShape3D;
	typedef stir::Reconstruction<Image3DF> Reconstruction3DF;
	typedef stir::IterativeReconstruction<Image3DF> IterativeReconstruction3DF;
	typedef stir::GeneralisedObjectiveFunction<Image3DF> ObjectiveFunction3DF;
	typedef stir::PoissonLogLikelihoodWithLinearModelForMean < Image3DF >
		PoissonLogLhLinModMean3DF;
	//PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF>
	typedef stir::ProjectorByBinPairUsingProjMatrixByBin ProjectorPairUsingMatrix;
	typedef stir::ProjMatrixByBinUsingRayTracing RayTracingMatrix;
	typedef stir::GeneralisedPrior<Image3DF> Prior3DF;
	typedef stir::QuadraticPrior<float> QuadPrior3DF;
	typedef stir::DataProcessor<Image3DF> DataProcessor3DF;
	typedef stir::TruncateToCylindricalFOVImageProcessor<float> CylindricFilter3DF;

}

#endif