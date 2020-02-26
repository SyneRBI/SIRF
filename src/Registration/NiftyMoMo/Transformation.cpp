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
#include "sirf/NiftyMoMo/Transformation.h"
#include "_reg_resampling.h"
#include "_reg_maths.h"
#include <iostream>

using namespace NiftyMoMo;


//--------------------------------
// Transformation::TransformImage
//--------------------------------
nifti_image* Transformation::TransformImage( nifti_image* sourceImgIn, nifti_image* targetImgIn )
{
#ifndef NDEBUG
  std::cout << "Called Transformation::TransformImage()" << std::endl;
#endif
  nifti_image* warpedImage;

  /// \todo Warped image is created here for now, however, we need to think about ownership of this image and deletion on destruction. 
  ///       Proposal: Generate an image class that takes care of this when it is no longer needed.
  warpedImage = nifti_copy_nim_info( targetImgIn );
  warpedImage->data = malloc( warpedImage->nvox * warpedImage->nbyper );

  nifti_image* dvfImg = this->GetDeformationVectorField( targetImgIn );
  reg_resampleImage( sourceImgIn, 
                     warpedImage,
                     dvfImg, 
                     nullptr, 
                     this->interpolation, 
                     this->warpedPaddingValue );
  
  return warpedImage;
}




//---------------------------------------
// Transformation::TransformImageAdjoint
//---------------------------------------
void Transformation::TransformImageAdjoint( const nifti_image * const sourceImage, const nifti_image * const sourceWeightsImage,
  const nifti_image * warpedImage, const nifti_image * const warpedWeightsImage )
{
  // Need to calculate the DVF for each point of the source image
  // note the member variable will be used, so ignore return value
  this->GetDeformationVectorField( sourceImage );

  // need to determine the interpolation method and dimensionality to call correct function
  
  // cubic interpolation
  if (3 == this->interpolation)
  {
    // 3D
    if (3 == this->deformationVectorFieldImage->nu) 
    { 
      this->CubicSplineTransformImageAdjoint3D( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
    else
    {
      // 2D
      this->CubicSplineTransformImageAdjoint2D( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
  }
  // nearest neighbour interpolation
  else if (0 == this->interpolation) 
  {
    if (3 == this->deformationVectorFieldImage->nu)
    {
      // 3D
      this->NearestNeighbourTransformImageAdjoint3D( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
    else 
    {
      // 2D
      this->NearestNeighbourTransformImageAdjoint2D( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
  }
  // tri-/bi-linear interpolation by default
  else 
  {
    if (3 == this->deformationVectorFieldImage->nu)
    {
      // 3D
      this->TrilinearTransformImageAdjoint( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
    else
    {
      // 2D
      this->BilinearTransformImageAdjoint( sourceImage, sourceWeightsImage, warpedImage, warpedWeightsImage );
    }
  }
}



//----------------------------------------
// Transformation::GetImageGradientWRTDVF
//----------------------------------------
void Transformation::GetImageGradientWRTDVF( nifti_image* sourceImage, nifti_image* outWarpedGradientImage )
{
#ifndef NDEBUG
  std::cout << "Called Transformation::GetImageGradientWRTDVF()" << std::endl;
#endif

  if (this->deformationVectorFieldImage == nullptr)
  {
    this->GetDeformationVectorField(sourceImage);
  }

  reg_getImageGradient( sourceImage, 
                        outWarpedGradientImage, 
                        this->deformationVectorFieldImage, 
                        nullptr, 
                        this->interpolation, 
                        this->warpedPaddingValue, 0 );

  return;
}

void Transformation::setDVF( nifti_image * DVF )
{
    if (this->deformationVectorFieldImage != nullptr && this->needToDeleteDVF)
    {
        nifti_image_free( this->deformationVectorFieldImage );
    }
    this->deformationVectorFieldImage = DVF;
    this->dvfImageUpdateRequired = false;
    // no need to delete as the external source should take care of that
    this->needToDeleteDVF = false;
}



//---------------------------------------------
// Transformation::CheckDVFImageUpdateRequired
//---------------------------------------------
bool Transformation::CheckDVFImageUpdateRequired(const nifti_image * const targetImageIn ) const
{
#ifndef NDEBUG
  std::cout << "Called Transformation::CheckDVFImageUpdateRequired()" << std::endl;
#endif

  if ( true == this->dvfImageUpdateRequired )
  {
    return true;
  }

  if ( nullptr == this->deformationVectorFieldImage )
  {
    return true;
  }

  // Check the geometry of the DVF and the target image
  if ( targetImageIn->nx != this->deformationVectorFieldImage->nx ) return true;
  if ( targetImageIn->ny != this->deformationVectorFieldImage->ny ) return true;
  if ( targetImageIn->nz != this->deformationVectorFieldImage->nz ) return true;
  if ( targetImageIn->sform_code != this->deformationVectorFieldImage->sform_code ) return true;
  if ( targetImageIn->qform_code != this->deformationVectorFieldImage->qform_code ) return true;

  // Check sform-matrix or qform matrix otherwise
  if ( targetImageIn->sform_code > 0 )
  {
    // Assuming consistency between sto_ijk and sto_xyz
    if ( fabs(targetImageIn->sto_ijk.m[0][0] - this->deformationVectorFieldImage->sto_ijk.m[0][0]) > 1.e-4f) return true;
    if ( fabs(targetImageIn->sto_ijk.m[0][1] - this->deformationVectorFieldImage->sto_ijk.m[0][1]) > 1.e-4f) return true;
    if ( fabs(targetImageIn->sto_ijk.m[0][2] - this->deformationVectorFieldImage->sto_ijk.m[0][2]) > 1.e-4f) return true;
    if ( fabs(targetImageIn->sto_ijk.m[0][3] - this->deformationVectorFieldImage->sto_ijk.m[0][3]) > 1.e-4f) return true;

    if ( fabs(targetImageIn->sto_ijk.m[1][0] - this->deformationVectorFieldImage->sto_ijk.m[1][0]) > 1.e-4f) return true;
    if ( fabs(targetImageIn->sto_ijk.m[1][1] - this->deformationVectorFieldImage->sto_ijk.m[1][1]) > 1.e-4f) return true;
    if ( fabs(targetImageIn->sto_ijk.m[1][2] - this->deformationVectorFieldImage->sto_ijk.m[1][2]) > 1.e-4f) return true;
    if ( fabs(targetImageIn->sto_ijk.m[1][3] - this->deformationVectorFieldImage->sto_ijk.m[1][3]) > 1.e-4f) return true;

    if ( fabs(targetImageIn->sto_ijk.m[2][0] - this->deformationVectorFieldImage->sto_ijk.m[2][0]) > 1.e-4f) return true;
    if ( fabs(targetImageIn->sto_ijk.m[2][1] - this->deformationVectorFieldImage->sto_ijk.m[2][1]) > 1.e-4f) return true;
    if ( fabs(targetImageIn->sto_ijk.m[2][2] - this->deformationVectorFieldImage->sto_ijk.m[2][2]) > 1.e-4f) return true;
    if ( fabs(targetImageIn->sto_ijk.m[2][3] - this->deformationVectorFieldImage->sto_ijk.m[2][3]) > 1.e-4f) return true;
  }
  else
  {
    // Assuming consistency between qto_ijk and qto_xyz
    if ( fabs(targetImageIn->qto_ijk.m[0][0] - this->deformationVectorFieldImage->qto_ijk.m[0][0]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[0][1] - this->deformationVectorFieldImage->qto_ijk.m[0][1]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[0][2] - this->deformationVectorFieldImage->qto_ijk.m[0][2]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[0][3] - this->deformationVectorFieldImage->qto_ijk.m[0][3]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[1][0] - this->deformationVectorFieldImage->qto_ijk.m[1][0]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[1][1] - this->deformationVectorFieldImage->qto_ijk.m[1][1]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[1][2] - this->deformationVectorFieldImage->qto_ijk.m[1][2]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[1][3] - this->deformationVectorFieldImage->qto_ijk.m[1][3]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[2][0] - this->deformationVectorFieldImage->qto_ijk.m[2][0]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[2][1] - this->deformationVectorFieldImage->qto_ijk.m[2][1]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[2][2] - this->deformationVectorFieldImage->qto_ijk.m[2][2]) > 1.e-4f ) return true;
    if ( fabs(targetImageIn->qto_ijk.m[2][3] - this->deformationVectorFieldImage->qto_ijk.m[2][3]) > 1.e-4f ) return true;
  }
  
  return false;
}




//----------------------------------------------------
// Transformation::CubicSplineTransformImageAdjoint3D
//----------------------------------------------------
void Transformation::CubicSplineTransformImageAdjoint3D( 
  const nifti_image * const sourceImage,
  const nifti_image * const sourceWeightsImage,
  const nifti_image * warpedImage,
  const nifti_image * const warpedWeightsImage )
{
  //some useful pointers
  PrecisionType *sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType *warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny * sourceImage->nz;
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of dynamic image not warped image
  PrecisionType *deformationFieldPtrZ = &deformationFieldPtrY[sourceVoxelNumber];

  //matrix from real space to warped image space
  const mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0)
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else 
    warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType xBasis[4], yBasis[4], zBasis[4], relative, world[3], position[3], *warpedWeiZPtr;
  PrecisionType *warpedWeiXYZPtr, yzBasis, xyzBasis, intensity;
  int index, a, b, c, firstVox[3], Y, Z;
  PrecisionType *warpedImgZPtr, *warpedImgXYZPtr, warpedImgValue;

  /*#if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, xBasis, yBasis, zBasis, relative, world, position, firstVox, \
    a, b, c, Y, Z, warpedImgZPtr, warpedWeiZPtr, warpedWeiXYZPtr, warpedImgXYZPtr, \
    yzBasis, xyzBasis, intensity, warpedImgValue) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, maskPtr, \
    warpedIJKMatrix, warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP*/

  //as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; index++)
  {
    //only push value if in mask and not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      //find where the def field maps this dynamic voxel to in the warped image
      //first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      world[2] = deformationFieldPtrZ[index];
      
      //now use matrix to go from world space to warped image space (i.e. static image space)
      reg_mat44_mul( warpedIJKMatrix, world, position );
      
      //round down to voxel before def field position
      firstVox[0] = static_cast<int>(reg_floor( position[0] ));
      firstVox[1] = static_cast<int>(reg_floor( position[1] ));
      firstVox[2] = static_cast<int>(reg_floor( position[2] ));
      
      //get cubic spline basis values for each axis
      // basis values along the x axis
      relative = position[0] - (PrecisionType)firstVox[0];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, xBasis );
      
      // basis values along the y axis
      relative = position[1] - (PrecisionType)firstVox[1];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, yBasis );
      
      // basis values along the z axis
      relative = position[2] - (PrecisionType)firstVox[2];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, zBasis );
      
      //as using cubic spline need two voxels before and two after def field position
      //so shift firstVox one more in -ve direction
      --firstVox[0];
      --firstVox[1];
      --firstVox[2];

      //loop over 4 voxels on z axis
      for (c = 0; c < 4; c++)
      {
        Z = firstVox[2] + c;
        if (-1 < Z && Z < warpedImage->nz)
        {
          warpedImgZPtr = &warpedImgPtr[Z * warpedImage->nx * warpedImage->ny];
          warpedWeiZPtr = &warpedWeightsPtr[Z * warpedWeightsImage->nx * warpedWeightsImage->ny];

          //loop over 4 voxels on y axis
          for (b = 0; b < 4; b++)
          {
            Y = firstVox[1] + b;
            if (-1 < Y && Y < warpedImage->ny)
            {
              warpedImgXYZPtr = &warpedImgZPtr[Y * warpedImage->nx + firstVox[0]];
              warpedWeiXYZPtr = &warpedWeiZPtr[Y * warpedWeightsImage->nx + firstVox[0]];
              yzBasis = yBasis[b] * zBasis[c];

              //loop over 4 voxels on x axis
              for (a = 0; a < 4; a++)
              {
                if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedImage->nx)
                {
                  xyzBasis = yzBasis * xBasis[a];
                  //#pragma omp atomic
                  warpedWeiXYZPtr[a] += sourceWeightsPtr[index] * xyzBasis;
                  intensity = ((PrecisionType)sourcePtr[index]) * xyzBasis;
                  
                  /// \todo: Decide if original switch statement is necessary at all. 
                  //switch (sourceImage->datatype)
                  //{
                  //case NIFTI_TYPE_FLOAT32:
                  //case NIFTI_TYPE_FLOAT64:
                  warpedImgValue = (PrecisionType)intensity;
                  //  break;
                  //case NIFTI_TYPE_UINT8:
                  //case NIFTI_TYPE_UINT16:
                  //case NIFTI_TYPE_UINT32:
                  //  warpedImgValue = (PrecisionType)(intensity > 0 ? reg_round( intensity ) : 0);
                  //  break;
                  //default:
                  //  warpedImgValue = (PrecisionType)reg_round( intensity );
                  //}//switch (dynamicImage->datatype)
//#pragma omp atomic
                  warpedImgXYZPtr[a] += warpedImgValue;
                }//if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedDynImage->nx)
              }//for (a = 0; a < 4; a++)
            }//if (-1 < Y && Y < warpedDynImage->ny)
          }//for (b = 0; b < 4; b++)
        }//if (-1 < Z && Z < warpedDynImage->nz)
      }//for (c = 0; c < 4; c++)
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}




//----------------------------------------------------
// Transformation::CubicSplineTransformImageAdjoint2D
//----------------------------------------------------
void Transformation::CubicSplineTransformImageAdjoint2D( 
  const nifti_image * const sourceImage,
  const nifti_image * const sourceWeightsImage,
  const nifti_image * warpedImage,
  const nifti_image * const warpedWeightsImage )
{
  //some useful pointers
  PrecisionType *sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType *warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny;
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of dynamic image not warped image
 
  //matrix from real space to warped image space
  const mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0)
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else 
    warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType xBasis[4], yBasis[4], relative, world[2], position[2], *warpedWeiXYPtr;
  PrecisionType xyBasis, intensity;
  int index, a, b, firstVox[2], Y;
  PrecisionType *warpedImgXYPtr, warpedImgValue;

  /*#if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, xBasis, yBasis, relative, world, position, firstVox, \
    a, b, Y, warpedWeiXYPtr, warpedImgXYPtr, xyBasis, intensity, warpedImgValue) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, maskPtr, warpedIJKMatrix, \
    warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP*/

  //as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; index++)
  {
    //only push value if in mask and not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      //find where the def field maps this dynamic voxel to in the warped image
      //first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      
      //now use matrix to go from world space to warped image space (i.e. static image space)
      position[0] = world[0] * (*warpedIJKMatrix).m[0][0] + world[1] * (*warpedIJKMatrix).m[0][1] + (*warpedIJKMatrix).m[0][3];
      position[1] = world[0] * (*warpedIJKMatrix).m[1][0] + world[1] * (*warpedIJKMatrix).m[1][1] + (*warpedIJKMatrix).m[1][3];
      
      //round down to voxel before def field position
      firstVox[0] = static_cast<int>(reg_floor( position[0] ));
      firstVox[1] = static_cast<int>(reg_floor( position[1] ));
      
      //get cubic spline basis values for each axis
      // basis values along the x axis
      relative = position[0] - (PrecisionType)firstVox[0];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, xBasis );
      
      // basis values along the y axis
      relative = position[1] - (PrecisionType)firstVox[1];
      relative = relative > 0 ? relative : 0;
      interpolantCubicSpline<PrecisionType>( relative, yBasis );
      
      //as using cubic spline need two voxels before and two after def field position
      //so shift firstVox one more in -ve direction
      --firstVox[0];
      --firstVox[1];

      //loop over 4 voxels on y axis
      for (b = 0; b < 4; b++)
      {
        Y = firstVox[1] + b;
        if (-1 < Y && Y < warpedImage->ny)
        {
          warpedImgXYPtr = &warpedImgPtr[Y * warpedImage->nx + firstVox[0]];
          warpedWeiXYPtr = &warpedWeightsPtr[Y * warpedWeightsImage->nx + firstVox[0]];

          //loop over 4 voxels on x axis
          for (a = 0; a < 4; a++)
          {
            if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedImage->nx)
            {
              xyBasis = yBasis[b] * xBasis[a];
              //#pragma omp atomic
              warpedWeiXYPtr[a] += sourceWeightsPtr[index] * xyBasis;
              intensity = ((PrecisionType)sourcePtr[index]) * xyBasis;
              /// \todo: Decide if original switch statement is necessary at all. 
              //switch (sourceImage->datatype)
              //{
              //case NIFTI_TYPE_FLOAT32:
              //case NIFTI_TYPE_FLOAT64:
              warpedImgValue = (PrecisionType)intensity;
              //  break;
              //case NIFTI_TYPE_UINT8:
              //case NIFTI_TYPE_UINT16:
              //case NIFTI_TYPE_UINT32:
              //  warpedImgValue = (ImageTYPE)(intensity > 0 ? reg_round( intensity ) : 0);
              //  break;
              //default:
              //  warpedImgValue = (ImageTYPE)reg_round( intensity );
              //}//switch (dynamicImage->datatype)
//#pragma omp atomic
              warpedImgXYPtr[a] += warpedImgValue;
            }//if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedDynImage->nx)
          }//for (a = 0; a < 4; a++)
        }//if (-1 < Y && Y < warpedDynImage->ny)
      }//for (b = 0; b < 4; b++)
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}




//---------------------------------------------------------
// Transformation::NearestNeighbourTransformImageAdjoint3D
//---------------------------------------------------------
void Transformation::NearestNeighbourTransformImageAdjoint3D( 
  const nifti_image * const sourceImage,
  const nifti_image * const sourceWeightsImage,
  const nifti_image * warpedImage,
  const nifti_image * const warpedWeightsImage )
{
  //some useful pointers
  PrecisionType *sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType *warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType *sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType *warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny * sourceImage->nz;
  PrecisionType *deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType *deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of source image not warped image
  PrecisionType *deformationFieldPtrZ = &deformationFieldPtrY[sourceVoxelNumber];

  //matrix from real space to warped image space
  const mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0)
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType world[3], position[3];
  int index, nearestVox[3], warpedInd;

  /*#if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, world, position, nearestVox, warpedInd) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, maskPtr, \
    warpedIJKMatrix, warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP*/

  //as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; index++)
  {
    //only push value if in mask and not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      //find where the def field maps this dynamic voxel to in the warped image
      //first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      world[2] = deformationFieldPtrZ[index];
      //now use matrix to go from world space to warped image space (i.e. static image space)
      reg_mat44_mul( warpedIJKMatrix, world, position );
      //round to nearest voxel to def field position
      nearestVox[0] = static_cast<int>(reg_round( position[0] ));
      nearestVox[1] = static_cast<int>(reg_round( position[1] ));
      nearestVox[2] = static_cast<int>(reg_round( position[2] ));

      //if nearestVox in image push value and weight
      if (-1 < nearestVox[0] && nearestVox[0] < warpedImage->nx &&
          -1 < nearestVox[1] && nearestVox[1] < warpedImage->ny &&
          -1 < nearestVox[2] && nearestVox[2] < warpedImage->nz)
      {
        warpedInd = (nearestVox[2] * warpedImage->ny + nearestVox[1]) * warpedImage->nx + nearestVox[0];
        //#pragma omp atomic
        warpedImgPtr[warpedInd] += sourcePtr[index];
        //#pragma omp atomic
        warpedWeightsPtr[warpedInd] += sourceWeightsPtr[index];
      }
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}




//---------------------------------------------------------
// Transformation::NearestNeighbourTransformImageAdjoint2D
//---------------------------------------------------------
void Transformation::NearestNeighbourTransformImageAdjoint2D( 
  const nifti_image * const sourceImage,
  const nifti_image * sourceWeightsImage,
  const nifti_image * const warpedImage,
  const nifti_image * const warpedWeightsImage )
{
  //some useful pointer definitions
  PrecisionType* sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType* warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType* sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType* warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny;
  PrecisionType* deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType* deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of dynamic image not warped image

  //matrix from real space to warped image space
  const mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0)
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else 
    warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType world[2], position[2];
  int index, nearestVox[2], warpedInd;

  /*#if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, world, position, nearestVox, warpedInd) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, maskPtr, \
    warpedIJKMatrix, warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP*/

  //as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; index++)
  {
    //only push value if in mask and not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      //find where the def field maps this dynamic voxel to in the warped image
      //first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      
      //now use matrix to go from world space to warped image space (i.e. static image space)
      position[0] = world[0] * (*warpedIJKMatrix).m[0][0] + world[1] * (*warpedIJKMatrix).m[0][1] + (*warpedIJKMatrix).m[0][3];
      position[1] = world[0] * (*warpedIJKMatrix).m[1][0] + world[1] * (*warpedIJKMatrix).m[1][1] + (*warpedIJKMatrix).m[1][3];
      
      //round to nearest voxel to def field position
      nearestVox[0] = static_cast<int>(reg_round( position[0] ));
      nearestVox[1] = static_cast<int>(reg_round( position[1] ));

      //if nearestVox in image push value and weight
      if (-1 < nearestVox[0] && nearestVox[0] < warpedImage->nx &&
          -1 < nearestVox[1] && nearestVox[1] < warpedImage->ny)
      {
        warpedInd = nearestVox[1] * warpedImage->nx + nearestVox[0];
        //#pragma omp atomic
        warpedImgPtr[warpedInd] += sourcePtr[index];
        //#pragma omp atomic
        warpedWeightsPtr[warpedInd] += sourceWeightsPtr[index];
      }
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}




//------------------------------------------------
// Transformation::TrilinearTransformImageAdjoint
//------------------------------------------------
void Transformation::TrilinearTransformImageAdjoint( 
  const nifti_image * const sourceImage,
  const nifti_image * const sourceWeightsImage,
  const nifti_image* warpedImage,
  const nifti_image * const warpedWeightsImage )
{
  //some useful pointer definitions
  PrecisionType* sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType* warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType* sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType* warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny * sourceImage->nz;
  PrecisionType* deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType* deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of source image not warped image
  PrecisionType* deformationFieldPtrZ = &deformationFieldPtrY[sourceVoxelNumber];
  
  //matrix from real space to warped image space
  const mat44 * warpedIJKMatrix;
  if (warpedImage->sform_code > 0) 
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else 
    warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType xBasis[2], yBasis[2], zBasis[2], relative, world[3], position[3], *warpedWeiZPtr;
  PrecisionType *warpedWeiXYZPtr, yzBasis, xyzBasis, intensity;
  int index, a, b, c, firstVox[3], Y, Z;
  PrecisionType *warpedImgZPtr, *warpedImgXYZPtr, warpedImgValue;

  /*#if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, xBasis, yBasis, zBasis, relative, world, position, firstVox, \
    a, b, c, Y, Z, warpedImgZPtr, warpedWeiZPtr, warpedWeiXYZPtr, warpedImgXYZPtr, \
    yzBasis, xyzBasis, intensity, warpedImgValue) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, deformationFieldPtrZ, maskPtr, \
    warpedIJKMatrix, warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP*/

  // as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; ++index)
  {
    // only push value if not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      // find where the def field maps this source voxel to in the warped image
      // first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      world[2] = deformationFieldPtrZ[index];

      // now use matrix to go from world space to warped image space (i.e. static image space)
      reg_mat44_mul( warpedIJKMatrix, world, position );
      
      // round down to voxel before def field position
      firstVox[0] = static_cast<int>(reg_floor( position[0] ));
      firstVox[1] = static_cast<int>(reg_floor( position[1] ));
      firstVox[2] = static_cast<int>(reg_floor( position[2] ));
      
      // get basis values for each axis
      // basis values along the x axis
      relative = position[0] - (PrecisionType)firstVox[0];
      xBasis[0] = (PrecisionType)(1.0 - relative);
      xBasis[1] = relative;
      // basis values along the y axis
      relative = position[1] - (PrecisionType)firstVox[1];
      yBasis[0] = (PrecisionType)(1.0 - relative);
      yBasis[1] = relative;
      // basis values along the z axis
      relative = position[2] - (PrecisionType)firstVox[2];
      zBasis[0] = (PrecisionType)(1.0 - relative);
      zBasis[1] = relative;

      //loop over 2 voxels on z axis
      for (c = 0; c < 2; c++)
      {
        Z = firstVox[2] + c;
        if (-1 < Z && Z < warpedImage->nz)
        {
          warpedImgZPtr = &warpedImgPtr[Z * warpedImage->nx * warpedImage->ny];
          warpedWeiZPtr = &warpedWeightsPtr[Z * warpedWeightsImage->nx * warpedWeightsImage->ny];

          //loop over 2 voxels on y axis
          for (b = 0; b < 2; b++)
          {
            Y = firstVox[1] + b;
            if (-1 < Y && Y < warpedImage->ny)
            {
              warpedImgXYZPtr = &warpedImgZPtr[Y * warpedImage->nx + firstVox[0]];
              warpedWeiXYZPtr = &warpedWeiZPtr[Y * warpedWeightsImage->nx + firstVox[0]];
              yzBasis = yBasis[b] * zBasis[c];

              //loop over 2 voxels on x axis
              for (a = 0; a < 2; a++)
              {
                if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedImage->nx)
                {
                  xyzBasis = yzBasis * xBasis[a];
                  
                  //#pragma omp atomic
                  warpedWeiXYZPtr[a] += sourceWeightsPtr[index] * xyzBasis;
                  intensity = ((PrecisionType)sourcePtr[index]) * xyzBasis;
                  /// \todo: Decide if original switch statement is necessary at all. 
                  //switch (sourceImage->datatype)
                  //{
                  //case NIFTI_TYPE_FLOAT32:
                  //case NIFTI_TYPE_FLOAT64:
                  warpedImgValue = (PrecisionType)intensity;
                  //  break;
                  //case NIFTI_TYPE_UINT8:
                  //case NIFTI_TYPE_UINT16:
                  //case NIFTI_TYPE_UINT32:
                  //  warpedImgValue = (PrecisionType)(intensity > 0 ? reg_round( intensity ) : 0);
                  //  break;
                  //default:
                  //  warpedImgValue = (PrecisionType)reg_round( intensity );
                  //}//switch (dynamicImage->datatype)
//#pragma omp atomic
                  warpedImgXYZPtr[a] += warpedImgValue;
                }//if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedDynImage->nx)
              }//for (a = 0; a < 2; a++)
            }//if (-1 < Y && Y < warpedDynImage->ny)
          }//for (b = 0; b < 2; b++)
        }//if (-1 < Z && Z < warpedDynImage->nz)
      }//for (c = 0; c < 2; c++)
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}




//-----------------------------------------------
// Transformation::BilinearTransformImageAdjoint
//-----------------------------------------------
void Transformation::BilinearTransformImageAdjoint( 
  const nifti_image* const sourceImage,
  const nifti_image* const sourceWeightsImage,
  const nifti_image* warpedImage,
  const nifti_image* const warpedWeightsImage )
{
  //some useful pointers
  PrecisionType* sourcePtr = static_cast<PrecisionType *>(sourceImage->data);
  PrecisionType* warpedImgPtr = static_cast<PrecisionType *>(warpedImage->data);
  PrecisionType* sourceWeightsPtr = static_cast<PrecisionType *>(sourceWeightsImage->data);
  PrecisionType* warpedWeightsPtr = static_cast<PrecisionType *>(warpedWeightsImage->data);
  size_t sourceVoxelNumber = (size_t)sourceImage->nx * sourceImage->ny;
  PrecisionType* deformationFieldPtrX = static_cast<PrecisionType *>(this->deformationVectorFieldImage->data);
  PrecisionType* deformationFieldPtrY = &deformationFieldPtrX[sourceVoxelNumber]; //as using push interpolation def field in space of dynamic image not warped image

  //matrix from real space to warped image space
  const mat44 *warpedIJKMatrix;
  if (warpedImage->sform_code > 0)
    warpedIJKMatrix = &(warpedImage->sto_ijk);
  else 
    warpedIJKMatrix = &(warpedImage->qto_ijk);

  //declare variables used in omp loop
  PrecisionType xBasis[2], yBasis[2], relative, world[2], position[2], *warpedWeiXYPtr;
  PrecisionType xyBasis, intensity;
  int index, a, b, firstVox[2], Y;
  PrecisionType *warpedImgXYPtr, warpedDynValue;

  /*#if defined (NDEBUG) && defined (_OPENMP)
  #pragma omp parallel for default(none) \
    private(index, xBasis, yBasis, relative, world, position, firstVox, \
    a, b, Y, warpedWeiXYPtr, warpedImgXYPtr, xyBasis, intensity, warpedImgValue) \
    shared(sourcePtr, warpedImgPtr, sourceWeightsPtr, warpedWeightsPtr, sourceVoxelNumber, \
    deformationFieldPtrX, deformationFieldPtrY, maskPtr, warpedIJKMatrix, \
    warpedDynImage, warpedWeiImage, dynamicImage)
  #endif // _OPENMP*/

  //as using push interpolation need to loop over voxels in dynamic image not warped image
  for (index = 0; index < sourceVoxelNumber; index++)
  {
    //only push value if in mask and not NaN
    if (sourcePtr[index] == sourcePtr[index])
    {
      //find where the def field maps this dynamic voxel to in the warped image
      //first find co-ord in world space
      world[0] = deformationFieldPtrX[index];
      world[1] = deformationFieldPtrY[index];
      //now use matrix to go from world space to warped image space (i.e. static image space)
      position[0] = world[0] * (*warpedIJKMatrix).m[0][0] + world[1] * (*warpedIJKMatrix).m[0][1] +
        (*warpedIJKMatrix).m[0][3];
      position[1] = world[0] * (*warpedIJKMatrix).m[1][0] + world[1] * (*warpedIJKMatrix).m[1][1] +
        (*warpedIJKMatrix).m[1][3];
      //round down to voxel before def field position
      firstVox[0] = static_cast<int>(reg_floor( position[0] ));
      firstVox[1] = static_cast<int>(reg_floor( position[1] ));
      //get cubic spline basis values for each axis
      // basis values along the x axis
      relative = position[0] - (PrecisionType)firstVox[0];
      xBasis[0] = (PrecisionType)(1.0 - relative);
      xBasis[1] = relative;
      // basis values along the y axis
      relative = position[1] - (PrecisionType)firstVox[1];
      yBasis[0] = (PrecisionType)(1.0 - relative);
      yBasis[1] = relative;

      //loop over 2 voxels on y axis
      for (b = 0; b < 2; b++)
      {
        Y = firstVox[1] + b;
        if (-1 < Y && Y < warpedImage->ny)
        {
          warpedImgXYPtr = &warpedImgPtr[Y*warpedImage->nx + firstVox[0]];
          warpedWeiXYPtr = &warpedWeightsPtr[Y*warpedWeightsImage->nx + firstVox[0]];

          //loop over 2 voxels on x axis
          for (a = 0; a < 2; a++)
          {
            if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedImage->nx)
            {
              xyBasis = yBasis[b] * xBasis[a];
              //#pragma omp atomic
              warpedWeiXYPtr[a] += sourceWeightsPtr[index] * xyBasis;
              intensity = ((PrecisionType)sourcePtr[index]) * xyBasis;
              /// \todo: Decide if original switch statement is necessary at all. 
              /// \todo: Remove reg_round macro with inline function
              //switch (sourceImage->datatype)
              //{
              //case NIFTI_TYPE_FLOAT32:
              //case NIFTI_TYPE_FLOAT64:
              warpedDynValue = (PrecisionType)intensity;
              //  break;
              //case NIFTI_TYPE_UINT8:
              //case NIFTI_TYPE_UINT16:
              //case NIFTI_TYPE_UINT32:
              //  warpedImgValue = (PrecisionType)(intensity > 0 ? reg_round( intensity ) : 0);
              //  break;
              //default:
              //  warpedImgValue = (PrecisionType)reg_round( intensity );
              //}//switch (dynamicImage->datatype)
//#pragma omp atomic
              warpedImgXYPtr[a] += warpedDynValue;
            }//if (-1 < (firstVox[0] + a) && (firstVox[0] + a) < warpedDynImage->nx)
          }//for (a = 0; a < 2; a++)
        }//if (-1 < Y && Y < warpedDynImage->ny)
      }//for (b = 0; b < 2; b++)
    }//if ((maskPtr[index])>-1 && sourcePtr[index] == sourcePtr[index])
  }//for (index = 0; index<sourceVoxelNumber; index++)
}



//----------------------------------------
// Transformation::ReorientateVectorImage
//----------------------------------------
void Transformation::ReorientateVectorImage(nifti_image* vectorFieldImageToReorientate, mat44 reorientationMatrix)
{
	PrecisionType *vecPtrX = static_cast<PrecisionType *>(vectorFieldImageToReorientate->data);
	PrecisionType *vecPtrY = &vecPtrX[vectorFieldImageToReorientate->nx*vectorFieldImageToReorientate->ny*vectorFieldImageToReorientate->nz*vectorFieldImageToReorientate->nt];
	PrecisionType *vecPtrZ = NULL;
	if (vectorFieldImageToReorientate->nu == 3)
		vecPtrZ = &vecPtrY[vectorFieldImageToReorientate->nx*vectorFieldImageToReorientate->ny*vectorFieldImageToReorientate->nz*vectorFieldImageToReorientate->nt];

	for (int n = 0; n < vectorFieldImageToReorientate->nx*vectorFieldImageToReorientate->ny*vectorFieldImageToReorientate->nz*vectorFieldImageToReorientate->nt; n++)
	{
		PrecisionType reorientatedValues[3] = { 0.0, 0.0, 0.0 };
		if (vecPtrZ == NULL) // 2D
		{
			reorientatedValues[0] =
				reorientationMatrix.m[0][0] * vecPtrX[n] +
				reorientationMatrix.m[1][0] * vecPtrY[n];
			reorientatedValues[1] =
				reorientationMatrix.m[0][1] * vecPtrX[n] +
				reorientationMatrix.m[1][1] * vecPtrY[n];
			vecPtrX[n] = reorientatedValues[0];
			vecPtrY[n] = reorientatedValues[1];
		}
		else // 3D
		{
			reorientatedValues[0] =
				reorientationMatrix.m[0][0] * vecPtrX[n] +
				reorientationMatrix.m[1][0] * vecPtrY[n] +
				reorientationMatrix.m[2][0] * vecPtrZ[n];
			reorientatedValues[1] =
				reorientationMatrix.m[0][1] * vecPtrX[n] +
				reorientationMatrix.m[1][1] * vecPtrY[n] +
				reorientationMatrix.m[2][1] * vecPtrZ[n];
			reorientatedValues[2] =
				reorientationMatrix.m[0][2] * vecPtrX[n] +
				reorientationMatrix.m[1][2] * vecPtrY[n] +
				reorientationMatrix.m[2][2] * vecPtrZ[n];
			vecPtrX[n] = reorientatedValues[0];
			vecPtrY[n] = reorientatedValues[1];
			vecPtrZ[n] = reorientatedValues[2];
		}
	}

	return;
}
