// ====================================================================================================   
//                                                                                                        
//   niftyMoMo: An implementation of the generalised motion modelling and image registration framework    
//                                                                                                        
//   Copyright (c) University College London (UCL). All rights reserved.                                  
//                                                                                                        
//   This software is distributed WITHOUT ANY WARRANTY; without even                                      
//   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR                                  
//   PURPOSE.                                                                                             
//                                                                                                        
//   See LICENSE.txt in the top level directory for details.                                              
//                                                                                                        
// ====================================================================================================   




//----------
// Includes 
//----------
#include "sirf/NiftyMoMo/BSplineTransformation.h"
#include "_reg_localTrans_regul.h"
#include "_reg_localTrans.h"

using namespace NiftyMoMo;

inline void nmm_exit(int Val, const char* fileName, int lineNumber)
{
  fprintf(stderr, "[niftyMoMo] Exit here. File: %s:%i\n", fileName, lineNumber);
  exit(Val);
};

inline void nmm_print_debug(const char* message){ printf("[niftyMoMo DEBUG] Message: %s\n", message); };
inline void nmm_print_warning(const char* message){ printf("[niftyMoMo WARNING] Message: %s\n", message); };
inline void nmm_print_error(const char* message){ fprintf(stderr, "[niftyMoMo ERROR] Message: %s\n", message); };

//----------------------------------------------
// BSplineTransformation::BSplineTransformation
//----------------------------------------------
BSplineTransformation::BSplineTransformation( const BSplineTransformation& transformToCopy )
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::BSplineTransformation( const BSplineTransformation& transformToCopy )" << std::endl;
#endif
  // Copy all parameters explicitly
  // B-spline member variables
  this->bendingEnergyWeight = transformToCopy.bendingEnergyWeight;
  this->linearEnergyWeight = transformToCopy.linearEnergyWeight;
  this->lastInitialisedLevel = transformToCopy.lastInitialisedLevel;
  this->warpedPaddingValue = transformToCopy.warpedPaddingValue;
  this->dvfImageUpdateRequired = transformToCopy.dvfImageUpdateRequired;

  // Copy the CPG-image
  this->controlPointGridImage = nifti_copy_nim_info( transformToCopy.controlPointGridImage );
  this->controlPointGridImage->data = malloc( controlPointGridImage->nvox * controlPointGridImage->nbyper );
  memcpy( this->controlPointGridImage->data, transformToCopy.controlPointGridImage->data, this->controlPointGridImage->nvox * this->controlPointGridImage->nbyper );
  
  // And DVF if it exists
  if ( transformToCopy.deformationVectorFieldImage != nullptr )
  {
    this->needToDeleteDVF = true; // always true since we're creating a copy
    this->deformationVectorFieldImage = nifti_copy_nim_info( transformToCopy.deformationVectorFieldImage );

    this->deformationVectorFieldImage->data = malloc( deformationVectorFieldImage->nvox * deformationVectorFieldImage->nbyper );
    memcpy( this->deformationVectorFieldImage->data, transformToCopy.deformationVectorFieldImage->data, this->deformationVectorFieldImage->nvox * this->deformationVectorFieldImage->nbyper );
  }
  else
  {
    this->needToDeleteDVF = false;
    this->deformationVectorFieldImage = nullptr;
  }

  // Base-class parameters
  this->numberOfParameters = transformToCopy.numberOfParameters;
  this->numberOfLevelsToPerform = transformToCopy.numberOfLevelsToPerform;
  this->interpolation = transformToCopy.interpolation;
}




//----------------------------------------------
// BSplineTransformation::BSplineTransformation
//----------------------------------------------
BSplineTransformation::BSplineTransformation( nifti_image* referenceImageIn,
                                              unsigned int numberOfLevelsToPerformIn,
                                              float* finalControlPointGridSpacing )
                                              :
                                              controlPointGridImage( nullptr ),
                                              bendingEnergyWeight( 0.f ),
                                              linearEnergyWeight( 0.f ),
                                              lastInitialisedLevel( -1 )
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::BSplineTransformation()" << std::endl;
#endif

  // Initialisation of baseclass parameters
  this->interpolation = 1;
  this->deformationVectorFieldImage = nullptr;
  this->warpedPaddingValue = std::numeric_limits<Transformation::PrecisionType>::quiet_NaN();
  this->dvfImageUpdateRequired = true;
  this->needToDeleteDVF = true;

  // Array needs to be initialised here for compatibility with VS2013
  // otherwise bSplineSpacing{0,0,0} in list above would have worked.
  float cpgSpacingInMM[3];
  cpgSpacingInMM[0] = finalControlPointGridSpacing[0];
  cpgSpacingInMM[1] = finalControlPointGridSpacing[1];
  cpgSpacingInMM[2] = finalControlPointGridSpacing[2];

  this->numberOfLevelsToPerform = numberOfLevelsToPerformIn;

  /* Convert the spacing from voxel to mm if necessary */
  //float spacingInMillimeter[3] = { finalControlPointGridSpacing[0], finalControlPointGridSpacing[1], finalControlPointGridSpacing[2] };
  if ( cpgSpacingInMM[0] < 0 ) cpgSpacingInMM[0] *= -1.0f * referenceImageIn->dx;
  if ( cpgSpacingInMM[1] < 0 ) cpgSpacingInMM[1] *= -1.0f * referenceImageIn->dy;
  if ( cpgSpacingInMM[2] < 0 ) cpgSpacingInMM[2] *= -1.0f * referenceImageIn->dz;

  // Define the spacing for the first level
  cpgSpacingInMM[0] *= powf( 2.0f, (float) (this->numberOfLevelsToPerform - 1) );
  cpgSpacingInMM[1] *= powf( 2.0f, (float) (this->numberOfLevelsToPerform - 1) );

  if ( referenceImageIn->nz > 1 )
    cpgSpacingInMM[2] *= powf( 2.0f, (float) (this->numberOfLevelsToPerform - 1) );
  else
    cpgSpacingInMM[2] *= 1.0f;

  this->numberOfLevelsToPerform = numberOfLevelsToPerformIn;

  // Create the control-point grid image
  reg_createControlPointGrid<PrecisionType>( &this->controlPointGridImage, referenceImageIn, cpgSpacingInMM );
  this->numberOfParameters = this->controlPointGridImage->nvox;

  return;
}




//-----------------------------------------------
// BSplineTransformation::~BSplineTransformation
//-----------------------------------------------
BSplineTransformation::~BSplineTransformation()
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::~BSplineTransformation()" << std::endl;
#endif

  // Clean the control-point grid image
  if (this->controlPointGridImage != nullptr)
  {
    nifti_image_free( this->controlPointGridImage );
    this->controlPointGridImage = nullptr;;
  }

  // Clean DVF image
  if ( this->deformationVectorFieldImage != nullptr && this->needToDeleteDVF)
  {
    nifti_image_free( this->deformationVectorFieldImage );
    this->deformationVectorFieldImage = nullptr;;
  }

  return;
}




//--------------------------------------
// BSplineTransformation::GetParameters
//--------------------------------------
BSplineTransformation::PrecisionType* BSplineTransformation::GetParameters()
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::GetParameters()" << std::endl;
#endif

  return (BSplineTransformation::PrecisionType*) this->controlPointGridImage->data;
}




//----------------------------------------------
// BSplineTransformation::GetNumberOfParameters
//----------------------------------------------
unsigned int BSplineTransformation::GetNumberOfParameters()
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::GetNumberOfParameters()" << std::endl;
#endif

  return static_cast<unsigned int>(this->controlPointGridImage->nvox);
}




//----------------------------------------
// BSplineTransformation::InitialiseLevel
//----------------------------------------
void BSplineTransformation::InitialiseLevel( unsigned int levelIn )
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::InitialiseLevel()" << std::endl;
#endif
  
  // Calculate the step size between the last initialised level and the requested one.
  const int stepBetweenLevels = static_cast<int>(levelIn) - this->lastInitialisedLevel;

  // Warn if level was previously initialised and don't do anything
  if (stepBetweenLevels == 0)
  {
    nmm_print_warning( "BSplineTransformation::InitialiseLevel()" );
    nmm_print_warning( "Level was already initialised. Neither refining grid, nor updating bending-energy weight." );
    return;
  }

  // Check that there is no jump larger than one between resolution levels
  if (stepBetweenLevels != 1)
  {
    nmm_print_error( "Initialisation only implemented to go from one level to the next." );
    nmm_exit( EXIT_FAILURE, __FILE__, __LINE__ );
  }

  // Only do a grid refinement on the subsequent levels
  if (0 != levelIn) reg_spline_refineControlPointGrid( this->controlPointGridImage );

  // update the weights of the constraint/regularisation term
  if (0 == levelIn)
  {
    this->bendingEnergyWeight = this->bendingEnergyWeight / powf( 16.f, (float)(this->numberOfLevelsToPerform) - 1.f );
  }
  else
  {
    this->bendingEnergyWeight *= 16.f;
  }
  
  // Remember which level was initialised last
  this->lastInitialisedLevel = levelIn;

  // Update the number of parameters
  this->numberOfParameters = this->controlPointGridImage->nvox;
}



//--------------------------------------------------
// BSplineTransformation::GetDeformationVectorField
//--------------------------------------------------
nifti_image* BSplineTransformation::GetDeformationVectorField( const nifti_image * const targetImageIn )
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::GetDeformationVectorField()" << std::endl;
#endif

  // Return the existing DVF it does not need to be updated
  if ( false == this->CheckDVFImageUpdateRequired( targetImageIn ) )
  {
    return this->deformationVectorFieldImage;
  }

  // Clear DVF image if it exists
  if ( nullptr != this->deformationVectorFieldImage && this->needToDeleteDVF )
  {
    nifti_image_free( this->deformationVectorFieldImage );
    this->deformationVectorFieldImage = nullptr;;
  }

  // Since we're creating our own DVF, we'll need to remember to delete it
  this->needToDeleteDVF = true;

  // Allocate the DVF
  // and fill it with the correct dimensions
  this->deformationVectorFieldImage = nifti_copy_nim_info( targetImageIn );
  this->deformationVectorFieldImage->dim[0] = this->deformationVectorFieldImage->ndim = 5;
  this->deformationVectorFieldImage->dim[4] = this->deformationVectorFieldImage->nt = 1;
  this->deformationVectorFieldImage->pixdim[4] = this->deformationVectorFieldImage->dt = 1.0;
  if ( this->controlPointGridImage->nu == 2 ) 
    this->deformationVectorFieldImage->dim[5] = this->deformationVectorFieldImage->nu = 2;
  else this->deformationVectorFieldImage->dim[5] = this->deformationVectorFieldImage->nu = 3;
  this->deformationVectorFieldImage->pixdim[5] = this->deformationVectorFieldImage->du = 1.0;
  this->deformationVectorFieldImage->dim[6] = this->deformationVectorFieldImage->nv = 1;
  this->deformationVectorFieldImage->pixdim[6] = this->deformationVectorFieldImage->dv = 1.0;
  this->deformationVectorFieldImage->dim[7] = this->deformationVectorFieldImage->nw = 1;
  this->deformationVectorFieldImage->pixdim[7] = this->deformationVectorFieldImage->dw = 1.0;
  this->deformationVectorFieldImage->nvox =
    (size_t)this->deformationVectorFieldImage->nx *
    (size_t)this->deformationVectorFieldImage->ny *
    (size_t)this->deformationVectorFieldImage->nz *
    (size_t)this->deformationVectorFieldImage->nt *
    (size_t)this->deformationVectorFieldImage->nu;
  this->deformationVectorFieldImage->nbyper = sizeof( PrecisionType );
  if ( sizeof( PrecisionType ) == sizeof( float ) )
    this->deformationVectorFieldImage->datatype = NIFTI_TYPE_FLOAT32;
  else this->deformationVectorFieldImage->datatype = NIFTI_TYPE_FLOAT64;
  this->deformationVectorFieldImage->data = (void *) calloc( this->deformationVectorFieldImage->nvox,
                                                             this->deformationVectorFieldImage->nbyper );
  this->deformationVectorFieldImage->intent_code = NIFTI_INTENT_VECTOR;
  memset( this->deformationVectorFieldImage->intent_name, 0, 16 );
  strcpy( this->deformationVectorFieldImage->intent_name, "NREG_TRANS" );
  this->deformationVectorFieldImage->intent_p1 = DEF_FIELD;
  this->deformationVectorFieldImage->scl_slope = 1.f;
  this->deformationVectorFieldImage->scl_inter = 0.f;

  if ( this->controlPointGridImage->num_ext > 0 )
  {
    nmm_print_error("BSplineTransformation::GetDeformationVectorField()");
    nmm_print_error( "The spline control point grid includes affine transformations." );
    nmm_print_error( "This function cannot be used when the control point grid includes affine transformations" );
    nmm_print_error( "as it assumes the axes of the control point grid and dense vector field are aligned." );
    nmm_exit( 1, __FILE__, __LINE__ );
  }

  //3D or 2D version?
  //checks that CPG and DVF are both 3D or both 2D
  if ( this->controlPointGridImage->nu == 2 && this->deformationVectorFieldImage->nu == 2 )
  {
    this->CPGToDVF2D( this->deformationVectorFieldImage );
   }
  else if ( this->controlPointGridImage->nu == 3 && this->deformationVectorFieldImage->nu == 3 )
  {
    this->CPGToDVF3D( this->deformationVectorFieldImage );
  }
  else
  {
    nmm_print_error( "reg_spline_CPGToDVF()" );
    nmm_print_error( "The spline control point grid and the dense vector field must have the same number of components for each vector" );
    nmm_print_error( "and these must be either 2 (for a 2D transformation) or 3 (for a 3D transformation)" );
    nmm_exit( 1, __FILE__, __LINE__ );
  }
  
  // Once the DVF image is calculated, remember that re-calculation only needs to happen if
  // 1) The parameters of the transformation are changed or 
  // 2) The target image has changed
  this->dvfImageUpdateRequired = false;
  
  return this->deformationVectorFieldImage;
}




//-----------------------------------
// BSplineTransformation::CPGToDVF2D
//-----------------------------------
void BSplineTransformation::CPGToDVF2D( nifti_image *denseVectorFieldImage )
{
  // pointers to values in the CPG and DVF image
  PrecisionType *CPGPtrX = static_cast<PrecisionType *>(this->controlPointGridImage->data);
  PrecisionType *CPGPtrY = &CPGPtrX[this->controlPointGridImage->nx * this->controlPointGridImage->ny];

  PrecisionType *DVFPtrX = static_cast<PrecisionType *>(denseVectorFieldImage->data);
  PrecisionType *DVFPtrY = &DVFPtrX[denseVectorFieldImage->nx * denseVectorFieldImage->ny];

  //matrix from DVF voxels to CPG voxels
  mat44 DVFToCPGMatrix, DVFXYZMatrix, CPGIJKMatrix;
  if ( denseVectorFieldImage->sform_code > 0 )
    DVFXYZMatrix = denseVectorFieldImage->sto_xyz;
  else
    DVFXYZMatrix = denseVectorFieldImage->qto_xyz;
  if ( this->controlPointGridImage->sform_code > 0 )
    CPGIJKMatrix = this->controlPointGridImage->sto_ijk;
  else
    CPGIJKMatrix = this->controlPointGridImage->qto_ijk;
  DVFToCPGMatrix = reg_mat44_mul( &(CPGIJKMatrix), &(DVFXYZMatrix) );
  PrecisionType DVFSpacingInCPG[2], DVFOriginInCPG[2];
  DVFSpacingInCPG[0] = DVFToCPGMatrix.m[0][0];
  DVFSpacingInCPG[1] = DVFToCPGMatrix.m[1][1];
  DVFOriginInCPG[0] = DVFToCPGMatrix.m[0][3];
  DVFOriginInCPG[1] = DVFToCPGMatrix.m[1][3];

  //declare variables used in omp loops
  int yDVF, yCPG, xDVF, xCPG, firstYCPG, firstXCPG, indCPG, indDVF;
  PrecisionType tmpBasis, bSplineBasis[4], *tmpXValsX, *tmpXValsY, valX, valY;


  //use openmp if available and not debug build
#if defined (NDEBUG) && defined (_OPENMP) && _OPENMP>201012
#pragma omp parallel for default(none) \
	shared(denseVectorFieldImage, DVFSpacingInCPG, \
	DVFOriginInCPG, CPGPtrX, CPGPtrY, DVFPtrX, DVFPtrY) \
	private(yDVF, yCPG, xDVF, xCPG, tmpBasis, bSplineBasis, tmpXValsX, tmpXValsY, \
	firstYCPG, firstXCPG, indCPG, indDVF, valX, valY)
#endif

  //loop over y in DVF
  for ( yDVF = 0; yDVF < denseVectorFieldImage->ny; yDVF++ )
  {
    //calc first index into DVF for this value of yDVF
    indDVF = yDVF * denseVectorFieldImage->nx;

    //find first yCPG to affect this yDVF
    tmpBasis = (PrecisionType) yDVF * DVFSpacingInCPG[1] + DVFOriginInCPG[1];
    firstYCPG = (int) tmpBasis - 1;

    //check if yCPGs used for this yDVF are in CPG
    if ( firstYCPG < 0 || firstYCPG + 4 > this->controlPointGridImage->ny )
    {
      //if outside CPG then assign NaN to all dense vector field values for this yDVF
      //loop over x in DVF
      for ( xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++ )
      {
        DVFPtrX[indDVF] = std::numeric_limits<PrecisionType>::quiet_NaN();
        DVFPtrY[indDVF] = std::numeric_limits<PrecisionType>::quiet_NaN();
        indDVF++;
      }
    }//one or more yCPGs used for this yDVF are outside CPG
    else
    {
      //all yCPGs for this yDVF are in CPG

      //find b-spline basis values for this yDVF
      tmpBasis = tmpBasis - (PrecisionType) firstYCPG - 1;
      get_BSplineBasisValues<PrecisionType>( tmpBasis, bSplineBasis );

      //initialise variables to store intermediate results
      tmpXValsX = new PrecisionType[this->controlPointGridImage->nx];
      tmpXValsY = new PrecisionType[this->controlPointGridImage->nx];

      //loop over x in CPG
      for ( xCPG = 0; xCPG < this->controlPointGridImage->nx; xCPG++ )
      {
        //sum the 4 CPs in y direction that affect yDVF
        valX = valY = 0.0;
        for ( yCPG = 0; yCPG < 4; yCPG++ )
        {
          indCPG = xCPG + (yCPG + firstYCPG) * this->controlPointGridImage->nx;
          valX += CPGPtrX[indCPG] * bSplineBasis[yCPG];
          valY += CPGPtrY[indCPG] * bSplineBasis[yCPG];
        }
        tmpXValsX[xCPG] = valX;
        tmpXValsY[xCPG] = valY;
      }// loop over x in CPG

      //loop over x in DVF
      for ( xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++ )
      {
        valX = valY = 0.0;

        //find first xCPG to affect this xDVF
        tmpBasis = (PrecisionType) xDVF * DVFSpacingInCPG[0] + DVFOriginInCPG[0];
        firstXCPG = (int) tmpBasis - 1;

        //find b-spline basis values for this xDVF
        tmpBasis = tmpBasis - (PrecisionType) firstXCPG - 1;
        get_BSplineBasisValues<PrecisionType>( tmpBasis, bSplineBasis );

        //check if xCPGs in CPG
        if ( firstXCPG < 0 || firstXCPG + 4 > this->controlPointGridImage->nx )
        {
          //if outside CPG then assign NaN to dense vector field
          valX = std::numeric_limits<PrecisionType>::quiet_NaN();
          valY = std::numeric_limits<PrecisionType>::quiet_NaN();
        }
        else
        {
          //sum the 4 CPs in the x direction that affect xDVF
          for ( xCPG = 0; xCPG < 4; xCPG++ )
          {
            valX += tmpXValsX[firstXCPG + xCPG] * bSplineBasis[xCPG];
            valY += tmpXValsY[firstXCPG + xCPG] * bSplineBasis[xCPG];
          }
        }

        //assigne DVF values
        DVFPtrX[indDVF] = valX;
        DVFPtrY[indDVF] = valY;

        //increment dvf index
        indDVF++;
      }//loop over x in DVF

      //delete intermediate results
      delete[] tmpXValsX;
      delete[] tmpXValsY;
    }//yCPGs used for this yDVF are all inside CPG
  }//loop over y in DVF
  
  return;
}




//-----------------------------------
// BSplineTransformation::CPGToDVF3D
//-----------------------------------
void BSplineTransformation::CPGToDVF3D( nifti_image *denseVectorFieldImage )
{
  //some useful pointers
  PrecisionType *CPGPtrX = static_cast<PrecisionType *>(this->controlPointGridImage->data);
  PrecisionType *CPGPtrY = &CPGPtrX[this->controlPointGridImage->nx * this->controlPointGridImage->ny * this->controlPointGridImage->nz];
  PrecisionType *CPGPtrZ = &CPGPtrY[this->controlPointGridImage->nx * this->controlPointGridImage->ny * this->controlPointGridImage->nz];

  PrecisionType *DVFPtrX = static_cast<PrecisionType *>(denseVectorFieldImage->data);
  PrecisionType *DVFPtrY = &DVFPtrX[denseVectorFieldImage->nx * denseVectorFieldImage->ny * denseVectorFieldImage->nz];
  PrecisionType *DVFPtrZ = &DVFPtrY[denseVectorFieldImage->nx * denseVectorFieldImage->ny * denseVectorFieldImage->nz];

  //matrix from DVF voxels to CPG voxels
  mat44 DVFToCPGMatrix, DVFXYZMatrix, CPGIJKMatrix;
  if ( denseVectorFieldImage->sform_code > 0 )
    DVFXYZMatrix = denseVectorFieldImage->sto_xyz;
  else
    DVFXYZMatrix = denseVectorFieldImage->qto_xyz;
  if ( this->controlPointGridImage->sform_code > 0 )
    CPGIJKMatrix = this->controlPointGridImage->sto_ijk;
  else
    CPGIJKMatrix = this->controlPointGridImage->qto_ijk;
  DVFToCPGMatrix = reg_mat44_mul( &(CPGIJKMatrix), &(DVFXYZMatrix) );
  PrecisionType DVFSpacingInCPG[3], DVFOriginInCPG[3];
  DVFSpacingInCPG[0] = DVFToCPGMatrix.m[0][0];
  DVFSpacingInCPG[1] = DVFToCPGMatrix.m[1][1];
  DVFSpacingInCPG[2] = DVFToCPGMatrix.m[2][2];
  DVFOriginInCPG[0] = DVFToCPGMatrix.m[0][3];
  DVFOriginInCPG[1] = DVFToCPGMatrix.m[1][3];
  DVFOriginInCPG[2] = DVFToCPGMatrix.m[2][3];

  //declare variables used in omp loops
  int zDVF, zCPG, yDVF, yCPG, xDVF, xCPG, firstZCPG, firstYCPG, firstXCPG, indCPG, indDVF, indTmp;
  PrecisionType tmpBasis, bSplineBasis[4], *tmpXYValsX, *tmpXYValsY, *tmpXYValsZ, *tmpXValsX, *tmpXValsY, *tmpXValsZ, valX, valY, valZ;


  //use openmp if available and not debug build
#if defined (NDEBUG) && defined (_OPENMP) && _OPENMP>201012
#pragma omp parallel for default(none) \
	shared(denseVectorFieldImage, DVFSpacingInCPG, \
	DVFOriginInCPG, CPGPtrX, CPGPtrY, CPGPtrZ, DVFPtrX, DVFPtrY, DVFPtrZ) \
	private(zDVF, zCPG, yDVF, yCPG, xDVF, xCPG, tmpBasis, bSplineBasis, \
	tmpXYValsX, tmpXYValsY, tmpXYValsZ, tmpXValsX, tmpXValsY, tmpXValsZ, \
	firstZCPG, firstYCPG, firstXCPG, indCPG, indDVF, indTmp, valX, valY, valZ)
#endif

  //loop over z in DVF
  for ( zDVF = 0; zDVF < denseVectorFieldImage->nz; zDVF++ )
  {
    //calc first index into DVF for this value of zDVF
    indDVF = zDVF * denseVectorFieldImage->nx * denseVectorFieldImage->ny;

    //find first zCPG to affect this zDVF
    tmpBasis = (PrecisionType) zDVF * DVFSpacingInCPG[2] + DVFOriginInCPG[2];
    firstZCPG = (int) tmpBasis - 1;

    //check if zCPGs used for this zDVF are in CPG
    if ( firstZCPG < 0 || firstZCPG + 4 > this->controlPointGridImage->nz )
    {
      //if outside CPG then assign NaN to all dense vector field values for this zDVF
      //loop over y and x in DVF
      for ( yDVF = 0; yDVF < denseVectorFieldImage->ny; yDVF++ )
      {
        for ( xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++ )
        {
          DVFPtrX[indDVF] = std::numeric_limits<PrecisionType>::quiet_NaN();
          DVFPtrY[indDVF] = std::numeric_limits<PrecisionType>::quiet_NaN();
          DVFPtrZ[indDVF] = std::numeric_limits<PrecisionType>::quiet_NaN();
          indDVF++;
        }
      }
    }//one or more zCPGs used for this zDVF are outside CPG
    else
    {
      //all zCPGs for this zDVF are in CPG

      //find b-spline basis values for this zDVF
      tmpBasis = tmpBasis - (PrecisionType) firstZCPG - 1;
      get_BSplineBasisValues<PrecisionType>( tmpBasis, bSplineBasis );

      //initialise variables to store intermediate results
      tmpXYValsX = new PrecisionType[this->controlPointGridImage->nx * this->controlPointGridImage->ny];
      tmpXYValsY = new PrecisionType[this->controlPointGridImage->nx * this->controlPointGridImage->ny];
      tmpXYValsZ = new PrecisionType[this->controlPointGridImage->nx * this->controlPointGridImage->ny];

      //loop over y in CPG
      for ( yCPG = 0; yCPG < this->controlPointGridImage->ny; yCPG++ )
      {
        //loop over x in CPG
        for ( xCPG = 0; xCPG < this->controlPointGridImage->nx; xCPG++ )
        {
          valX = valY = valZ = 0.0;
          indTmp = xCPG + yCPG * this->controlPointGridImage->nx;
          //sum the 4 CPs in z direction that affect zDVF
          for ( zCPG = 0; zCPG < 4; zCPG++ )
          {
            indCPG = indTmp + (zCPG + firstZCPG) * this->controlPointGridImage->nx * this->controlPointGridImage->ny;
            valX += CPGPtrX[indCPG] * bSplineBasis[zCPG];
            valY += CPGPtrY[indCPG] * bSplineBasis[zCPG];
            valZ += CPGPtrZ[indCPG] * bSplineBasis[zCPG];
          }
          tmpXYValsX[indTmp] = valX;
          tmpXYValsY[indTmp] = valY;
          tmpXYValsZ[indTmp] = valZ;
        }//x in CPG
      }//y in CPG

      //loop over y in DVF
      for ( yDVF = 0; yDVF < denseVectorFieldImage->ny; yDVF++ )
      {
        //find first yCPG to affect this yDVF
        tmpBasis = (PrecisionType) yDVF * DVFSpacingInCPG[1] + DVFOriginInCPG[1];
        firstYCPG = (int) tmpBasis - 1;

        //check if yCPGs used for this yDVF are in CPG
        if ( firstYCPG < 0 || firstYCPG + 4 > this->controlPointGridImage->ny )
        {
          //if outside CPG then assign NaN to all dense vector field values for this yDVF

          //calc index into DVF for this value of yDVF and zDVF
          indDVF = (zDVF * denseVectorFieldImage->ny + yDVF) * denseVectorFieldImage->nx;
          //loop over x in DVF
          for ( xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++ )
          {
            DVFPtrX[indDVF] = std::numeric_limits<PrecisionType>::quiet_NaN();
            DVFPtrY[indDVF] = std::numeric_limits<PrecisionType>::quiet_NaN();
            DVFPtrZ[indDVF] = std::numeric_limits<PrecisionType>::quiet_NaN();
            indDVF++;
          }
        }//one or more yCPGs used for this yDVF are outside CPG
        else
        {
          //all yCPGs for this yDVF are in CPG

          //find b-spline basis values for this yDVF
          tmpBasis = tmpBasis - (PrecisionType) firstYCPG - 1;
          get_BSplineBasisValues<PrecisionType>( tmpBasis, bSplineBasis );

          //initialise variables to store intermediate results
          tmpXValsX = new PrecisionType[this->controlPointGridImage->nx];
          tmpXValsY = new PrecisionType[this->controlPointGridImage->nx];
          tmpXValsZ = new PrecisionType[this->controlPointGridImage->nx];

          //loop over x in CPG
          for ( xCPG = 0; xCPG < this->controlPointGridImage->nx; xCPG++ )
          {
            //sum the 4 CPs in y direction that affect yDVF
            valX = valY = valZ = 0.0;
            for ( yCPG = 0; yCPG < 4; yCPG++ )
            {
              indTmp = xCPG + (yCPG + firstYCPG) * this->controlPointGridImage->nx;
              valX += tmpXYValsX[indTmp] * bSplineBasis[yCPG];
              valY += tmpXYValsY[indTmp] * bSplineBasis[yCPG];
              valZ += tmpXYValsZ[indTmp] * bSplineBasis[yCPG];
            }
            tmpXValsX[xCPG] = valX;
            tmpXValsY[xCPG] = valY;
            tmpXValsZ[xCPG] = valZ;
          }// loop over x in CPG

          //loop over x in DVF
          for ( xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++ )
          {
            valX = valY = valZ = 0.0;

            //find first xCPG to affect this xDVF
            tmpBasis = (PrecisionType) xDVF * DVFSpacingInCPG[0] + DVFOriginInCPG[0];
            firstXCPG = (int) tmpBasis - 1;

            //find b-spline basis values for this xDVF
            tmpBasis = tmpBasis - (PrecisionType) firstXCPG - 1;
            get_BSplineBasisValues<PrecisionType>( tmpBasis, bSplineBasis );

            //check if xCPGs in CPG
            if ( firstXCPG < 0 || firstXCPG + 4 > this->controlPointGridImage->nx )
            {
              //if outside CPG then assign NaN to dense vector field
              valX = std::numeric_limits<PrecisionType>::quiet_NaN();
              valY = std::numeric_limits<PrecisionType>::quiet_NaN();
              valZ = std::numeric_limits<PrecisionType>::quiet_NaN();
            }
            else
            {
              //sum the 4 CPs in the x direction that affect xDVF
              for ( xCPG = 0; xCPG < 4; xCPG++ )
              {
                valX += tmpXValsX[firstXCPG + xCPG] * bSplineBasis[xCPG];
                valY += tmpXValsY[firstXCPG + xCPG] * bSplineBasis[xCPG];
                valZ += tmpXValsZ[firstXCPG + xCPG] * bSplineBasis[xCPG];
              }
            }

            //assigne DVF values
            DVFPtrX[indDVF] = valX;
            DVFPtrY[indDVF] = valY;
            DVFPtrZ[indDVF] = valZ;

            //increment dvf index
            indDVF++;
          }//loop over x in DVF

          //delete intermediate results
          delete[] tmpXValsX;
          delete[] tmpXValsY;
          delete[] tmpXValsZ;

        }//yCPGs used for this yDVF are all inside CPG
      }//loop over y in DVF

      //delete intermediate results
      delete[] tmpXYValsX;
      delete[] tmpXYValsY;
      delete[] tmpXYValsZ;

    }//zCPGs used for this zDVF are all inside CPG
  }//loop over z in DVF

  return;
}




//-----------------------------------
// BSplineTransformation::DVFToCPG2D
//-----------------------------------
void BSplineTransformation::DVFToCPG2D( nifti_image *controlPointGridImage,
                                        nifti_image *denseVectorFieldImage,
                                        bool normaliseKernel )
{
  //some useful pointers
  PrecisionType *CPGPtrX = static_cast<PrecisionType *>(controlPointGridImage->data);
  PrecisionType *CPGPtrY = &CPGPtrX[controlPointGridImage->nx * controlPointGridImage->ny];

  PrecisionType *DVFPtrX = static_cast<PrecisionType *>(denseVectorFieldImage->data);
  PrecisionType *DVFPtrY = &DVFPtrX[denseVectorFieldImage->nx * denseVectorFieldImage->ny];

  //matrix from CPG voxels to DVF voxels
  mat44 CPGToDVFMatrix, CPGXYZMatrix, DVFIJKMatrix;
  if (denseVectorFieldImage->sform_code > 0)
    DVFIJKMatrix = denseVectorFieldImage->sto_ijk;
  else
    DVFIJKMatrix = denseVectorFieldImage->qto_ijk;
  if (controlPointGridImage->sform_code > 0)
    CPGXYZMatrix = controlPointGridImage->sto_xyz;
  else
    CPGXYZMatrix = controlPointGridImage->qto_xyz;
  CPGToDVFMatrix = reg_mat44_mul( &(DVFIJKMatrix), &(CPGXYZMatrix) );
  PrecisionType CPGSpacingInDVF[2], CPGOriginInDVF[2], kernelWidth[2];
  CPGSpacingInDVF[0] = CPGToDVFMatrix.m[0][0];
  CPGSpacingInDVF[1] = CPGToDVFMatrix.m[1][1];
  CPGOriginInDVF[0] = CPGToDVFMatrix.m[0][3];
  CPGOriginInDVF[1] = CPGToDVFMatrix.m[1][3];
  kernelWidth[0] = 2 * CPGSpacingInDVF[0];
  kernelWidth[1] = 2 * CPGSpacingInDVF[1];

  //initialise CPG to 0
  reg_tools_multiplyValueToImage( controlPointGridImage, controlPointGridImage, 0.0 );

  //declare variables used in omp loops
  int yDVF, yCPG, xDVF, xCPG, firstYDVF, lastYDVF, firstXDVF, lastXDVF, indCPG, indDVF;
  PrecisionType yCPGinDVF, xCPGinDVF, bSplineBasis, relative, kernelSum, *tmpXValsX, *tmpXValsY;


  //use openmp if available and not debug build
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
	shared(controlPointGridImage, denseVectorFieldImage, CPGSpacingInDVF, \
	CPGOriginInDVF, kernelWidth, CPGPtrX, CPGPtrY, DVFPtrX, DVFPtrY, normaliseKernel) \
	private(yDVF, yCPG, xDVF, xCPG, firstYDVF, lastYDVF, firstXDVF, lastXDVF, \
	indCPG, indDVF, yCPGinDVF, xCPGinDVF, bSplineBasis, relative, kernelSum, \
	tmpXValsX, tmpXValsY)
#endif

  //loop over y in CPG
  for (yCPG = 0; yCPG < controlPointGridImage->ny; yCPG++)
  {
    //find coord of yCPG in DVF, and first+last yDVF to affect this yCPG
    yCPGinDVF = (PrecisionType)yCPG * CPGSpacingInDVF[1] + CPGOriginInDVF[1];
    firstYDVF = (int)(ceil( yCPGinDVF - kernelWidth[1] ));
    lastYDVF = (int)(floor( yCPGinDVF + kernelWidth[1] ));

    //check if all yDVFs outside DVF
    if (lastYDVF < 0 || firstYDVF >= denseVectorFieldImage->ny)
    {
      //all yDVF outside DVF so set CPG = 0 for this yCPG

      //calc index into CPG for this yCPG
      indCPG = yCPG*controlPointGridImage->nx;

      //loop over x in CPG
      for (xCPG = 0; xCPG < controlPointGridImage->nx; xCPG++)
      {
        //set CPG to 0
        CPGPtrX[indCPG] = 0.0;
        CPGPtrY[indCPG++] = 0.0;
      }
    }//all yDVF outside DVF
    else
    {
      //check if some yDVFs outside DVF
      if (firstYDVF < 0)
        firstYDVF = 0;
      if (lastYDVF >= denseVectorFieldImage->ny)
        lastYDVF = denseVectorFieldImage->ny - 1;

      //initialise variables to store intermediate results
      tmpXValsX = (PrecisionType *)calloc( denseVectorFieldImage->nx, sizeof( PrecisionType ) );
      tmpXValsY = (PrecisionType *)calloc( denseVectorFieldImage->nx, sizeof( PrecisionType ) );
      kernelSum = 0.0;

      //calc index into DVF for firstYDVF
      indDVF = firstYDVF*denseVectorFieldImage->nx;

      //loop over yDVF that affect this yCPG
      for (yDVF = firstYDVF; yDVF <= lastYDVF; yDVF++)
      {
        //calculate b-spline basis
        relative = fabs( (PrecisionType)(yDVF - yCPGinDVF) / CPGSpacingInDVF[1] );
        if (relative < 1.0)
          bSplineBasis = (3 * relative*relative*(relative - 2) + 4.0f) / 6.0f;
        else if (relative < 2.0)
        {
          relative = 2.0f - relative;
          bSplineBasis = relative*relative*relative / 6.0f;
        }
        else
        {
          bSplineBasis = 0.0;
        }
        kernelSum += bSplineBasis;

        //loop over x in DVF
        for (xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++)
        {
          //store temporary result
          tmpXValsX[xDVF] += bSplineBasis*DVFPtrX[indDVF];
          tmpXValsY[xDVF] += bSplineBasis*DVFPtrY[indDVF++];
        }
      }//loop over yDVF that affect this yCPG

      //if normalise kernel divide temp results by kernelSum
      if (normaliseKernel)
      {
        //loop over x in DVF
        for (xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++)
        {
          //divide by kernel sum
          tmpXValsX[xDVF] /= kernelSum;
          tmpXValsY[xDVF] /= kernelSum;
        }
      }
      //calc index into CPG for this yCPG
      indCPG = yCPG*controlPointGridImage->nx;

      //loop over x in CPG
      for (xCPG = 0; xCPG < controlPointGridImage->nx; xCPG++)
      {
        //find coord of xCPG in DVF, and first+last xDVF to affect this xCPG
        xCPGinDVF = (PrecisionType)xCPG * CPGSpacingInDVF[0] + CPGOriginInDVF[0];
        firstXDVF = (int)(ceil( xCPGinDVF - kernelWidth[0] ));
        lastXDVF = (int)(floor( xCPGinDVF + kernelWidth[0] ));

        //check if all xDVFs outside DVF
        if (lastXDVF < 0 || firstXDVF >= denseVectorFieldImage->nx)
        {
          //all xDVF outside DVF so set CPG = 0 for this xCPG
          CPGPtrX[indCPG] = 0.0;
          CPGPtrY[indCPG++] = 0.0;
        }
        else
        {
          //check if some xDVFs outside DVF
          if (firstXDVF < 0)
            firstXDVF = 0;
          if (lastXDVF >= denseVectorFieldImage->nx)
            lastXDVF = denseVectorFieldImage->nx - 1;

          kernelSum = 0.0;

          //loop over xDVF that affect this xCPG
          for (xDVF = firstXDVF; xDVF <= lastXDVF; xDVF++)
          {
            //calculate b-spline basis
            relative = fabs( (PrecisionType)(xDVF - xCPGinDVF) / CPGSpacingInDVF[0] );
            if (relative < 1.0)
              bSplineBasis = (3 * relative*relative*(relative - 2) + 4.0f) / 6.0f;
            else if (relative < 2.0)
            {
              relative = 2.0f - relative;
              bSplineBasis = relative*relative*relative / 6.0f;
            }
            else
            {
              bSplineBasis = 0.0;
            }
            kernelSum += bSplineBasis;

            //store result
            CPGPtrX[indCPG] += bSplineBasis*tmpXValsX[xDVF];
            CPGPtrY[indCPG] += bSplineBasis*tmpXValsY[xDVF];
          }//loop over xDVF that affect this xCPG

          //if normalise kernel divide result by kernelSum
          if (normaliseKernel)
          {
            CPGPtrX[indCPG] /= kernelSum;
            CPGPtrY[indCPG] /= kernelSum;
          }

          //increment index into CPG
          indCPG++;
        }//check if all xDVFs outside DVF
      }//loop over x in CPG
      //free memory used for temporary X results
      free( tmpXValsX );
      free( tmpXValsY );
    }//check if all yDVFs outside DVF
  }//loop over y in CPG
}




//-----------------------------------
// BSplineTransformation::DVFToCPG2D
//-----------------------------------
void BSplineTransformation::DVFToCPG3D( nifti_image *controlPointGridImage,
                                        nifti_image *denseVectorFieldImage,
                                        bool normaliseKernel )
{
  //some useful pointers
  PrecisionType *CPGPtrX = static_cast<PrecisionType *>(controlPointGridImage->data);
  PrecisionType *CPGPtrY = &CPGPtrX[controlPointGridImage->nx * controlPointGridImage->ny * controlPointGridImage->nz];
  PrecisionType *CPGPtrZ = &CPGPtrY[controlPointGridImage->nx * controlPointGridImage->ny * controlPointGridImage->nz];

  PrecisionType *DVFPtrX = static_cast<PrecisionType *>(denseVectorFieldImage->data);
  PrecisionType *DVFPtrY = &DVFPtrX[denseVectorFieldImage->nx * denseVectorFieldImage->ny * denseVectorFieldImage->nz];
  PrecisionType *DVFPtrZ = &DVFPtrY[denseVectorFieldImage->nx * denseVectorFieldImage->ny * denseVectorFieldImage->nz];

  //matrix from CPG voxels to DVF voxels
  mat44 CPGToDVFMatrix, CPGXYZMatrix, DVFIJKMatrix;
  if (denseVectorFieldImage->sform_code > 0)
    DVFIJKMatrix = denseVectorFieldImage->sto_ijk;
  else
    DVFIJKMatrix = denseVectorFieldImage->qto_ijk;
  if (controlPointGridImage->sform_code > 0)
    CPGXYZMatrix = controlPointGridImage->sto_xyz;
  else
    CPGXYZMatrix = controlPointGridImage->qto_xyz;
  CPGToDVFMatrix = reg_mat44_mul( &(DVFIJKMatrix), &(CPGXYZMatrix) );
  PrecisionType CPGSpacingInDVF[3], CPGOriginInDVF[3], kernelWidth[3];
  CPGSpacingInDVF[0] = CPGToDVFMatrix.m[0][0];
  CPGSpacingInDVF[1] = CPGToDVFMatrix.m[1][1];
  CPGSpacingInDVF[2] = CPGToDVFMatrix.m[2][2];
  CPGOriginInDVF[0] = CPGToDVFMatrix.m[0][3];
  CPGOriginInDVF[1] = CPGToDVFMatrix.m[1][3];
  CPGOriginInDVF[2] = CPGToDVFMatrix.m[2][3];
  kernelWidth[0] = fabs( 2 * CPGSpacingInDVF[0] );
  kernelWidth[1] = fabs( 2 * CPGSpacingInDVF[1] );
  kernelWidth[2] = fabs( 2 * CPGSpacingInDVF[2] );

  //initialise CPG to 0
  reg_tools_multiplyValueToImage( controlPointGridImage, controlPointGridImage, 0.0 );

  //declare variables used in omp loops
  int zDVF, zCPG, yDVF, yCPG, xDVF, xCPG, firstZDVF, lastZDVF, firstYDVF, lastYDVF, \
    firstXDVF, lastXDVF, indCPG, indDVF, indTmp;
  PrecisionType zCPGinDVF, yCPGinDVF, xCPGinDVF, bSplineBasis, relative, kernelSum, \
    *tmpXYValsX, *tmpXYValsY, *tmpXYValsZ, *tmpXValsX, *tmpXValsY, *tmpXValsZ, \
    valX, valY, valZ;


  //use openmp if available and not debug build
#if defined (NDEBUG) && defined (_OPENMP)
#pragma omp parallel for default(none) \
	shared(controlPointGridImage, denseVectorFieldImage, CPGSpacingInDVF, \
	CPGOriginInDVF, kernelWidth, CPGPtrX, CPGPtrY, CPGPtrZ, DVFPtrX, DVFPtrY, \
	DVFPtrZ, normaliseKernel) \
	private(zDVF, zCPG, yDVF, yCPG, xDVF, xCPG, firstZDVF, lastZDVF, firstYDVF, \
	lastYDVF, firstXDVF, lastXDVF, indCPG, indDVF, zCPGinDVF, yCPGinDVF, xCPGinDVF,\
	bSplineBasis, relative, kernelSum, tmpXYValsX, tmpXYValsY, tmpXYValsZ, tmpXValsX, \
	tmpXValsY, tmpXValsZ, valX, valY, valZ, indTmp)
#endif

  //loop over z in CPG
  for (zCPG = 0; zCPG < controlPointGridImage->nz; zCPG++)
  {
    //find coord of zCPG in DVF, and first+last zDVF to affect this zCPG
    zCPGinDVF = (PrecisionType)zCPG * CPGSpacingInDVF[2] + CPGOriginInDVF[2];
    firstZDVF = (int)(ceil( zCPGinDVF - kernelWidth[2] ));
    lastZDVF = (int)(floor( zCPGinDVF + kernelWidth[2] ));

    //check if all zDVFs outside DVF
    if (lastZDVF < 0 || firstZDVF >= denseVectorFieldImage->nz)
    {
      //all zDVF outside DVF so set CPG = 0 for this zCPG

      //calc index into CPG for this zCPG
      indCPG = zCPG*controlPointGridImage->nx*controlPointGridImage->ny;

      //loop over all y and x in CPG
      for (yCPG = 0; yCPG < controlPointGridImage->ny; yCPG++)
      {
        for (xCPG = 0; xCPG < controlPointGridImage->nx; xCPG++)
        {
          CPGPtrX[indCPG] = 0.0;
          CPGPtrY[indCPG] = 0.0;
          CPGPtrZ[indCPG++] = 0.0;
        }
      }
    }//all zDVFs outside DVF
    else
    {
      //check if some zDVFs outside DVF
      if (firstZDVF < 0)
        firstZDVF = 0;
      if (lastZDVF >= denseVectorFieldImage->nz)
        lastZDVF = denseVectorFieldImage->nz - 1;

      //initialise variables to store intermediate results
      tmpXYValsX = (PrecisionType *)calloc( denseVectorFieldImage->nx * denseVectorFieldImage->ny, sizeof( PrecisionType ) );
      tmpXYValsY = (PrecisionType *)calloc( denseVectorFieldImage->nx * denseVectorFieldImage->ny, sizeof( PrecisionType ) );
      tmpXYValsZ = (PrecisionType *)calloc( denseVectorFieldImage->nx * denseVectorFieldImage->ny, sizeof( PrecisionType ) );
      kernelSum = 0.0;

      //calc index into DVF for firstZDVF
      indDVF = firstZDVF * denseVectorFieldImage->nx * denseVectorFieldImage->ny;

      //loop over zDVF that affect this zCPG
      for (zDVF = firstZDVF; zDVF <= lastZDVF; zDVF++)
      {
        //calculate b-spline basis
        relative = fabs( (PrecisionType)(zDVF - zCPGinDVF) / CPGSpacingInDVF[2] );
        if (relative < 1.0)
          bSplineBasis = (3 * relative*relative*(relative - 2) + 4.0) / 6.0;
        else if (relative < 2.0)
        {
          relative = 2.0 - relative;
          bSplineBasis = relative*relative*relative / 6.0;
        }
        else
        {
          bSplineBasis = 0.0;
        }
        kernelSum += bSplineBasis;

        //loop over y in DVF
        indTmp = 0;
        for (yDVF = 0; yDVF < denseVectorFieldImage->ny; yDVF++)
        {
          //loop over x in DVF
          for (xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++)
          {
            //store temporary result
            tmpXYValsX[indTmp] += bSplineBasis*DVFPtrX[indDVF];
            tmpXYValsY[indTmp] += bSplineBasis*DVFPtrY[indDVF];
            tmpXYValsZ[indTmp++] += bSplineBasis*DVFPtrZ[indDVF++];
          }
        }
      }//loop over zDVF that affect this zCPG

      //if normalise kernel divide temp results by kernelSum
      if (normaliseKernel)
      {
        //loop over y in DVF
        indTmp = 0;
        for (yDVF = 0; yDVF < denseVectorFieldImage->ny; yDVF++)
        {
          //loop over x in DVF
          for (xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++)
          {
            //divide by kernel sum
            tmpXYValsX[indTmp] /= kernelSum;
            tmpXYValsY[indTmp] /= kernelSum;
            tmpXYValsZ[indTmp++] /= kernelSum;
          }
        }
      }


      //calc index into CPG for this zCPG
      indCPG = zCPG*controlPointGridImage->nx*controlPointGridImage->ny;

      //loop over y in CPG
      for (yCPG = 0; yCPG < controlPointGridImage->ny; yCPG++)
      {
        //find coord of yCPG in DVF, and first+last yDVF to affect this yCPG
        yCPGinDVF = (PrecisionType)yCPG * CPGSpacingInDVF[1] + CPGOriginInDVF[1];
        firstYDVF = (int)(ceil( yCPGinDVF - kernelWidth[1] ));
        lastYDVF = (int)(floor( yCPGinDVF + kernelWidth[1] ));

        //check if all yDVFs outside DVF
        if (lastYDVF < 0 || firstYDVF >= denseVectorFieldImage->ny)
        {
          //all yDVFs outside DVF so set CPG = 0 for this yCPG

          //loop over all x in CPG
          for (xCPG = 0; xCPG < controlPointGridImage->nx; xCPG++)
          {
            CPGPtrX[indCPG] = 0.0;
            CPGPtrY[indCPG] = 0.0;
            CPGPtrZ[indCPG++] = 0.0;
          }
        }
        else
        {
          //check if some yDVFs outside DVF
          if (firstYDVF < 0)
            firstYDVF = 0;
          if (lastYDVF >= denseVectorFieldImage->ny)
            lastYDVF = denseVectorFieldImage->ny - 1;

          //initialise variables to store intermediate results
          tmpXValsX = (PrecisionType *)calloc( denseVectorFieldImage->nx, sizeof( PrecisionType ) );
          tmpXValsY = (PrecisionType *)calloc( denseVectorFieldImage->nx, sizeof( PrecisionType ) );
          tmpXValsZ = (PrecisionType *)calloc( denseVectorFieldImage->nx, sizeof( PrecisionType ) );
          kernelSum = 0.0;

          //calc index into temp results for first yDVF
          indTmp = firstYDVF*denseVectorFieldImage->nx;

          //loop over yDVF that affect this yCPG
          for (yDVF = firstYDVF; yDVF <= lastYDVF; yDVF++)
          {
            //calculate b-spline basis
            relative = fabs( (PrecisionType)(yDVF - yCPGinDVF) / CPGSpacingInDVF[1] );
            if (relative < 1.0)
              bSplineBasis = (3 * relative*relative*(relative - 2) + 4.0) / 6.0;
            else if (relative < 2.0)
            {
              relative = 2.0 - relative;
              bSplineBasis = relative*relative*relative / 6.0;
            }
            else
            {
              bSplineBasis = 0.0;
            }
            kernelSum += bSplineBasis;

            //loop over x in DVF
            for (xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++)
            {
              //store temporary result
              tmpXValsX[xDVF] += bSplineBasis*tmpXYValsX[indTmp];
              tmpXValsY[xDVF] += bSplineBasis*tmpXYValsY[indTmp];
              tmpXValsZ[xDVF] += bSplineBasis*tmpXYValsZ[indTmp++];
            }
          }//loop over yDVF that affect this yCPG

          //if normalise kernel divide temp results by kernelSum
          if (normaliseKernel)
          {
            //loop over x in DVF
            for (xDVF = 0; xDVF < denseVectorFieldImage->nx; xDVF++)
            {
              //divide by kernel sum
              tmpXValsX[xDVF] /= kernelSum;
              tmpXValsY[xDVF] /= kernelSum;
              tmpXValsZ[xDVF] /= kernelSum;
            }
          }

          //loop over x in CPG
          for (xCPG = 0; xCPG < controlPointGridImage->nx; xCPG++)
          {
            //find coord of xCPG in DVF, and first+last xDVF to affect this xCPG
            xCPGinDVF = (PrecisionType)xCPG * CPGSpacingInDVF[0] + CPGOriginInDVF[0];
            firstXDVF = (int)(ceil( xCPGinDVF - kernelWidth[0] ));
            lastXDVF = (int)(floor( xCPGinDVF + kernelWidth[0] ));

            //check if all xDVFs outside DVF
            if (lastXDVF < 0 || firstXDVF >= denseVectorFieldImage->nx)
            {
              CPGPtrX[indCPG] = 0.0;
              CPGPtrY[indCPG] = 0.0;
              CPGPtrZ[indCPG++] = 0.0;
            }
            else
            {
              //check if some xDVFs outside DVF
              if (firstXDVF < 0)
                firstXDVF = 0;
              if (lastXDVF >= denseVectorFieldImage->nx)
                lastXDVF = denseVectorFieldImage->nx - 1;

              kernelSum = 0.0;

              //loop over xDVF that affect this xCPG
              for (xDVF = firstXDVF; xDVF <= lastXDVF; xDVF++)
              {
                //calculate b-spline basis
                relative = fabs( (PrecisionType)(xDVF - xCPGinDVF) / CPGSpacingInDVF[0] );
                if (relative < 1.0)
                  bSplineBasis = (3 * relative*relative*(relative - 2) + 4.0) / 6.0;
                else if (relative < 2.0)
                {
                  relative = 2.0 - relative;
                  bSplineBasis = relative*relative*relative / 6.0;
                }
                else
                {
                  bSplineBasis = 0.0;
                }
                kernelSum += bSplineBasis;

                //store result
                CPGPtrX[indCPG] += bSplineBasis*tmpXValsX[xDVF];
                CPGPtrY[indCPG] += bSplineBasis*tmpXValsY[xDVF];
                CPGPtrZ[indCPG] += bSplineBasis*tmpXValsZ[xDVF];
              }//loop over xDVF that affect this xCPG

              //if normalise kernel divide result by kernelSum
              if (normaliseKernel)
              {
                CPGPtrX[indCPG] /= kernelSum;
                CPGPtrY[indCPG] /= kernelSum;
                CPGPtrZ[indCPG] /= kernelSum;
              }

              //increment index into CPG
              indCPG++;
            }//check if all zDVFs outside DVF
          }//loop over x in CPG
          //free memory used for temporary X results
          free( tmpXValsX );
          free( tmpXValsY );
          free( tmpXValsZ );
        }//check if all yDVFs outside DVF
      }//loop over y in CPG
      //free memory used for temporary XY results
      free( tmpXYValsX );
      free( tmpXYValsY );
      free( tmpXYValsZ );
    }//check if all zDVFs outside DVF
  }//loop over z in CPG
}


////-----------------------------------------------
//// BSplineTransformation::ReorientateVectorImage
////-----------------------------------------------
//void BSplineTransformation::ReorientateVectorImage( nifti_image* vectorFieldImageToReorientate, mat44 reorientationMatrix )
//{
//  PrecisionType *vecPtrX = static_cast<PrecisionType *>(vectorFieldImageToReorientate->data);
//  PrecisionType *vecPtrY = &vecPtrX[vectorFieldImageToReorientate->nx*vectorFieldImageToReorientate->ny*vectorFieldImageToReorientate->nz*vectorFieldImageToReorientate->nt];
//  PrecisionType *vecPtrZ = NULL;
//  if (vectorFieldImageToReorientate->nu == 3)
//    vecPtrZ = &vecPtrY[vectorFieldImageToReorientate->nx*vectorFieldImageToReorientate->ny*vectorFieldImageToReorientate->nz*vectorFieldImageToReorientate->nt];
//
//  for (int n = 0; n < vectorFieldImageToReorientate->nx*vectorFieldImageToReorientate->ny*vectorFieldImageToReorientate->nz*vectorFieldImageToReorientate->nt; n++)
//  {
//    PrecisionType reorientatedValues[3] = { 0.0, 0.0, 0.0 };
//    if (vecPtrZ == NULL) // 2D
//    {
//      reorientatedValues[0] =
//        reorientationMatrix.m[0][0] * vecPtrX[n] +
//        reorientationMatrix.m[1][0] * vecPtrY[n];
//      reorientatedValues[1] =
//        reorientationMatrix.m[0][1] * vecPtrX[n] +
//        reorientationMatrix.m[1][1] * vecPtrY[n];
//      vecPtrX[n] = reorientatedValues[0];
//      vecPtrY[n] = reorientatedValues[1];
//    }
//    else // 3D
//    {
//      reorientatedValues[0] =
//        reorientationMatrix.m[0][0] * vecPtrX[n] +
//        reorientationMatrix.m[1][0] * vecPtrY[n] +
//        reorientationMatrix.m[2][0] * vecPtrZ[n];
//      reorientatedValues[1] =
//        reorientationMatrix.m[0][1] * vecPtrX[n] +
//        reorientationMatrix.m[1][1] * vecPtrY[n] +
//        reorientationMatrix.m[2][1] * vecPtrZ[n];
//      reorientatedValues[2] =
//        reorientationMatrix.m[0][2] * vecPtrX[n] +
//        reorientationMatrix.m[1][2] * vecPtrY[n] +
//        reorientationMatrix.m[2][2] * vecPtrZ[n];
//      vecPtrX[n] = reorientatedValues[0];
//      vecPtrY[n] = reorientatedValues[1];
//      vecPtrZ[n] = reorientatedValues[2];
//    }
//  }
//  
//  return;
//}




//-------------------------------------------------------------------------
// BSplineTransformation::GetConstraintGradientWRTTransformationParameters
//-------------------------------------------------------------------------
BSplineTransformation::PrecisionType* BSplineTransformation::GetConstraintGradientWRTTransformationParameters()
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::GetRegularisationGradientWRTTransformationParameters()" << std::endl;
#endif
  
  // Generate the output image, same size as CPG
  nifti_image* constraintGradientImage = nifti_copy_nim_info( this->controlPointGridImage );
  constraintGradientImage->data = (PrecisionType*) calloc( constraintGradientImage->nvox, sizeof( PrecisionType ) );

  if (this->bendingEnergyWeight > 0)
  {
    reg_spline_approxBendingEnergyGradient( this->controlPointGridImage, constraintGradientImage, this->bendingEnergyWeight );
  }
  if (this->linearEnergyWeight > 0)
  {
    reg_spline_approxLinearEnergyGradient( this->controlPointGridImage, constraintGradientImage, this->linearEnergyWeight );
  }

  // Copy over the pointer to the data and delete the image
  PrecisionType* outConstraintGradientData = (PrecisionType*) constraintGradientImage->data;
  constraintGradientImage->data = nullptr;
  nifti_image_free( constraintGradientImage );

  return  outConstraintGradientData;
}




//-------------------------------------------
// BSplineTransformation::GetConstraintValue
//-------------------------------------------
double BSplineTransformation::GetConstraintValue()
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::GetConstraintValue()" << std::endl;
#endif

  double constraintSum=0.0;
  if ( this->bendingEnergyWeight > 0 )
  {
    constraintSum += this->bendingEnergyWeight * reg_spline_approxBendingEnergy( this->controlPointGridImage );
  }
  if ( this->linearEnergyWeight > 0 )
  {
    constraintSum += this->linearEnergyWeight * reg_spline_approxLinearEnergy( this->controlPointGridImage );
  }

  return constraintSum;
}




//------------------------------------------------------------------
// BSplineTransformation::GetDVFGradientWRTTransformationParameters
//------------------------------------------------------------------
BSplineTransformation::PrecisionType* BSplineTransformation::GetDVFGradientWRTTransformationParameters( nifti_image* denseDVFIn )
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::GetDVFGradientWRTTransformationParameters()" << std::endl;
#endif

  // Take the input DVF and project it to the transformation parametrs. 
  // This is essentially a convolution with the B-Spline basis function.

  // Generate the output data structure
  nifti_image* outDVFGradWRTTrafoParams = nifti_copy_nim_info( this->controlPointGridImage );
  outDVFGradWRTTrafoParams->data = malloc( outDVFGradWRTTrafoParams->nbyper * outDVFGradWRTTrafoParams->nvox );


  // 3D or 2D version?
  // check that CPG and DVF are both 3D or both 2D
  if (this->controlPointGridImage->nu == 2 && outDVFGradWRTTrafoParams->nu == 2)
  {
    this->DVFToCPG2D( outDVFGradWRTTrafoParams, denseDVFIn, false );
  }
  else if (this->controlPointGridImage->nu == 3 && outDVFGradWRTTrafoParams->nu == 3)
  {
    this->DVFToCPG3D( outDVFGradWRTTrafoParams, denseDVFIn, false );
  }
  else
  {
    nmm_print_error( "Error in BSplineTransformation::GetDVFGradientWRTTransformationParameters()" );
    nmm_print_error( "The spline control point grid and the dense vector field must have the same number of components for each vector" );
    nmm_print_error( "and these must be either 2 (for a 2D transformation) or 3 (for a 3D transformation)" );
    nmm_exit( 1, __FILE__, __LINE__ );
  }

  // Copy over the data pointer from the the image and detach it. Then delete the image.
  PrecisionType* outDVFGradWRTTrafoParamData = (PrecisionType*) outDVFGradWRTTrafoParams->data;
  outDVFGradWRTTrafoParams->data = nullptr;
  nifti_image_free( outDVFGradWRTTrafoParams );
  
  return outDVFGradWRTTrafoParamData;
}




//----------------------------------------------
// BSplineTransformation::SetLinearEnergyWeight
//----------------------------------------------
void BSplineTransformation::SetLinearEnergyWeight( double linearEnergyWeightIn )
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::SetLinearEnergyWeight()" << std::endl;
#endif

  if ( linearEnergyWeightIn < 0.f )
  {
    char msg[200];
    sprintf( msg, "Linear energy weight has to be positive." );
    nmm_print_error( msg );
    nmm_exit( 1, __FILE__, __LINE__ );
  }

  this->linearEnergyWeight = linearEnergyWeightIn;
  return;
}




//--------------------------------------
// BSplineTransformation::SetParameters
//--------------------------------------
void BSplineTransformation::SetParameters( BSplineTransformation::PrecisionType* parametersIn, bool parametersAreDisplacements )
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::SetParameters()" << std::endl;
#endif
  
  // Copy the parameters
  memcpy( this->controlPointGridImage->data, parametersIn, this->controlPointGridImage->nvox * this->controlPointGridImage->nbyper );
  
  // Make sure the dvf is updated if required
  this->dvfImageUpdateRequired = true;

  // Parameters were indicated to describe a displacement (relative) and NOT a deformation (absolute). 
  // Thus convert over so that DVF generation works correctly.
  if ( parametersAreDisplacements )
  {
    reg_getDeformationFromDisplacement( this->controlPointGridImage );
  }
}




//-----------------------------------------------
// BSplineTransformation::SetBendingEnergyWeight
//-----------------------------------------------
void BSplineTransformation::SetBendingEnergyWeight( double bendingEnergyWeightIn )
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::SetBendingEnergyWeight()" << std::endl;
#endif

  if (bendingEnergyWeightIn < 0.f)
  {
    char msg[200];
    sprintf( msg, "Linear energy weight has to be positive." );
    nmm_print_error( msg );
    nmm_exit( 1, __FILE__, __LINE__ );
  }

  this->bendingEnergyWeight = bendingEnergyWeightIn;
  return;
}




//-----------------------------------------------
// BSplineTransformation::SetBendingEnergyWeight
//-----------------------------------------------
BSplineTransformation::PrecisionType BSplineTransformation::GetSumOfPenaltyWeights()
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::GetSumOfPenaltyWeights()" << std::endl;
#endif
  
  // Note: The bending energy weight is updated from level to level. 
  //       However, when the weight is requested, the updating between levels is irrelevant.
  //       Thus we need to recover the original weight. 

  // If no level was initialised, then the 
  PrecisionType finalBendingEnergyWeight = this->bendingEnergyWeight;
  if (this->lastInitialisedLevel >= 0)
  {
    finalBendingEnergyWeight *= powf( 16.0, this->numberOfLevelsToPerform - this->lastInitialisedLevel - 1 );
  }
  PrecisionType sumOfPenaltyWeights = finalBendingEnergyWeight + this->linearEnergyWeight;
  
  return sumOfPenaltyWeights;
}



//------------------------------------------------------------
// BSplineTransformation::GetMaxTransformationParameterLength
//------------------------------------------------------------
BSplineTransformation::PrecisionType BSplineTransformation::GetMaxTransformationParameterLength( PrecisionType * parametersIn )
{
  // Calculate the total number of voxels in the control-point image
  unsigned int voxNumber = this->controlPointGridImage->nx * this->controlPointGridImage->ny * this->controlPointGridImage->nz;
  
  // Allocate the return value
   PrecisionType maxGradValue = (PrecisionType) 0;
  
  // Prepare pointers to the data
  PrecisionType* ptrX = parametersIn;
  PrecisionType *ptrY = &ptrX[voxNumber];
  
  // 3D case
  if (this->controlPointGridImage->nz > 1) {
    PrecisionType *ptrZ = &ptrY[voxNumber];
    for (unsigned int i = 0; i < voxNumber; ++i) {
      PrecisionType valX = 0, valY = 0, valZ = 0;
      valX = *ptrX++;
      valY = *ptrY++;
      valZ = *ptrZ++;
      PrecisionType length = (PrecisionType)(sqrt( valX * valX + valY * valY + valZ * valZ ));
      maxGradValue = (length > maxGradValue) ? length : maxGradValue;
    }
  } // otherwise assume 2D case
  else {
    for (int i = 0; i < voxNumber; i++) {
      PrecisionType valX = 0, valY = 0;
      valX = *ptrX++;
      valY = *ptrY++;
      PrecisionType length = (PrecisionType)(sqrt( valX * valX + valY * valY ));
      maxGradValue = (length > maxGradValue) ? length : maxGradValue;
    }
  }
  return maxGradValue;
}





//-------------------------------------------------
// BSplineTransformation::GetTransformationAsImage
//-------------------------------------------------
nifti_image * BSplineTransformation::GetTransformationAsImage()
{
	return this->controlPointGridImage;
}




//---------------------------------
// BSplineTransformation::DeepCopy
//---------------------------------
std::shared_ptr<Transformation> BSplineTransformation::DeepCopy()
{
#ifndef NDEBUG
  std::cout << "Called BSplineTransformation::DeepCopy()" << std::endl;
#endif

  std::shared_ptr<BSplineTransformation> copiedTransform = std::make_shared<BSplineTransformation>( *this );
  return copiedTransform;
}




//--------------------------------------------------------
// BSplineTransformation::DisplayTransformationParameters
//--------------------------------------------------------
void BSplineTransformation::DisplayTransformationParameters()
{
  printf( "BSpline-transformation parameters\n" );
  printf( "\t* control point grid image dimension: %i x %i x %i\n", 
          this->controlPointGridImage->nx,
          this->controlPointGridImage->ny,
          this->controlPointGridImage->nz );

  printf( "\t* image spacing: %g x %g x %g mm\n",
          this->controlPointGridImage->dx,
          this->controlPointGridImage->dy,
          this->controlPointGridImage->dz );

  printf( "\t* bending energy weight (at current level): %g \n", this->bendingEnergyWeight );
  printf( "\t* linear energy weight (at current level): %g \n", this->linearEnergyWeight );

#ifndef NDEBUG
  if (this->controlPointGridImage->sform_code > 0)
  {
    reg_mat44_disp( &(this->controlPointGridImage->sto_xyz), (char *)"\t* CPG sform matrix" );
  }
  else
  {
    reg_mat44_disp( &(this->controlPointGridImage->qto_xyz), (char *)"\t* CPG qform matrix" );
  }

#endif

  return;
}
