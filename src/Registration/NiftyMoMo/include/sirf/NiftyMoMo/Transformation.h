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



#pragma once


#include "nifti1_io.h"
#include <memory>

namespace NiftyMoMo {

template <class PrecisionType>
void interpolantCubicSpline( PrecisionType ratio, PrecisionType *basis )
{
  if (ratio < 0.0) ratio = 0.0; //reg_rounding error
  PrecisionType FF = ratio * ratio;
  basis[0] = (PrecisionType)((ratio * ((2.0 - ratio)*ratio - 1.0)) / 2.0);
  basis[1] = (PrecisionType)((FF * (3.0*ratio - 5.0) + 2.0) / 2.0);
  basis[2] = (PrecisionType)((ratio * ((4.0 - 3.0*ratio)*ratio + 1.0)) / 2.0);
  basis[3] = (PrecisionType)((ratio - 1.0) * FF / 2.0);
}

template <class PrecisionType>
void interpolantCubicSpline( PrecisionType ratio, PrecisionType *basis, PrecisionType *derivative )
{
  interpolantCubicSpline<PrecisionType>( ratio, basis );
  if (ratio < 0.0) ratio = 0.0; //reg_rounding error
  PrecisionType FF = ratio * ratio;
  derivative[0] = (PrecisionType)((4.0*ratio - 3.0*FF - 1.0) / 2.0);
  derivative[1] = (PrecisionType)((9.0*ratio - 10.0) * ratio / 2.0);
  derivative[2] = (PrecisionType)((8.0*ratio - 9.0*FF + 1) / 2.0);
  derivative[3] = (PrecisionType)((3.0*ratio - 2.0) * ratio / 2.0);
}


/** Transformation class
 *  Abstract base class that must be reimplemented to obtain a functioning transformation class.
 */
class Transformation
{
public:
  typedef float PrecisionType;

  virtual ~Transformation(){};

  /** Initialise the transformation and all internal parameters. 
   *  \param referenceImageIn The reference image defining the space where the transformation will be defined.
   */
  virtual void InitialiseLevel(unsigned int level) = 0;
  
  /** Get a pointer to the transformation parameters
   */
  virtual PrecisionType* GetParameters() = 0;

  /** Get the number of parameters required to describe this transformation completely.
   */
  virtual unsigned int GetNumberOfParameters() = 0;

  /** Set the transformation parameters
   */
  virtual void SetParameters( PrecisionType* paramsIn, bool parametersAreDisplacements ) = 0;

  /** Get the deformation vector field according to the current parameters
   *  \param targetImageIn Nifti image is used to allocate the DVF. The DVF will remain 
   *                       part of the transformation class and will be deleted if the 
   *                       transformation is deleted.
   */
  virtual nifti_image* GetDeformationVectorField( const nifti_image * const targetImageIn ) = 0;
  
  /** Compute the gradient of the DVF with respect to the transformation parameters
   *  \f$ \frac{\partial \textbf{DVF}_t}{\partial \mathbf{M_t}} \f$. 
   *  \param denseDVFIn pointer to the nifti image for which the best match of transformation 
   *                    parameters will be found
   *  \return pointer to the internally allocated output data of length numberOfParameters. Needs to be freed outside 
   *          of this class since no internal references will be kept.
   */
  virtual PrecisionType* GetDVFGradientWRTTransformationParameters( nifti_image* denseDVFIn ) = 0;
  
  /** Get the gradient of the regularisation/constraint with respect to the transformation parameters
   *  \f$ \frac{\partial \mathcal{R}_t}{\partial\mathbf{M}_t} \f$
   */
  virtual PrecisionType* GetConstraintGradientWRTTransformationParameters() = 0;
  
  /** Get the value of the constraint term
   */
  virtual double GetConstraintValue() = 0;
  
  /** Transform/warp an imgae
   */
  nifti_image* TransformImage( nifti_image* sourceImage, nifti_image* referenceImage );

  /** Transform an image using the "push-interpolation". This function was originally implemented
   *  specifically for reg_resp, hence the functionality will be provided by protected function of 
   *  this class.
   *  \param sourceImage Image that will be push-transformed into the warped image
   *  \param sourceWeightsImage Weights associated with the source image. This considers 
   *                            the image acquisition process where the acquired image 
   *                            may not have contributed equally to the voxel in the source 
   *                            image. 
   *  \param warpedImage Image into which the intensities will be pushed. Note that 
   *                     intensities will be accumulated.
   *  \param warpedWeightsImage Warped (and accumulated) weights of the source image.
   */
  void TransformImageAdjoint(const nifti_image * const sourceImage, const nifti_image * const sourceWeightsImage,
    const nifti_image *warpedImage, const nifti_image * const warpedWeightsImage );

  /** Generate a deep copy of the transformation and return it.
   */
  virtual std::shared_ptr<Transformation> DeepCopy() = 0;

  /** Print the current transformation parameters to the command line.
   */
  virtual void DisplayTransformationParameters() = 0;

  /** Calculate the sum of the penalty weights
   */
  virtual PrecisionType GetSumOfPenaltyWeights() = 0;

  /** Calculate the gradient of an image with respect to the deformation vector field.
   *  This calculates \f$ \frac{\partial \mathbf{I}_{T_t}}{\partial \textbf{DVF}_t} \f$ and is
   *  not dependent on any parametrisation of a derived transformation. 
   *  \param sourceImage Nifti-image structure with the image from which the derivative will be calculated.
   *  \param outWarpedGradientImage Nifti image structure wich will be filled with the calculated gradient. 
                                    Has to be allocated outside of this function.
   */
  virtual void GetImageGradientWRTDVF(nifti_image* sourceImage, nifti_image* outWarpedGradientImage);

  /** Function to reorientate vectors in 5D nifti image.
   *  Vectors which are stored along the 5th (u) dimension will be reoriented
   *  \param vectorFieldToReorientate Vector image to be reorientated
   *  \param reorientationMatrix 4 x 4 Matrix used to reorientate vectors
   */
  void ReorientateVectorImage(nifti_image* vectorFieldToReorientate, mat44 reorientationMatrix);


  /** Function to determine the maximum length of a given set of parameters. 
   *   This length depends on how the transformation parameters correspond internally to each other, i.e. how the x, y, and z
   *   component are stored.
   * 
   */
  virtual PrecisionType GetMaxTransformationParameterLength( PrecisionType* parametersIn ) = 0;


  /** Function that returns the current transformation as an image. This is not the DVF, but an image of the parameters instead. 
   *  Hence it can only be implemented by any derived classes. 
   */
  virtual nifti_image* GetTransformationAsImage() = 0;


  /** Get the padding value.
   */
  inline PrecisionType GetPaddingValue() { return this->warpedPaddingValue; };

  /** Set the padding value.
   */
  void SetPaddingValue(const PrecisionType paddingValue) { this->warpedPaddingValue = paddingValue; };


  /// Set interpolation
  void set_interpolation(const int interpolation_) { interpolation = interpolation_; }

  /// Set DVF
  void setDVF(nifti_image *DVF);

protected:

  /** Function to check if the DVF needs to be updated according to the target image
   *  \param targetImageIn The geometry if the DVF is checked against the target image
   */
  bool CheckDVFImageUpdateRequired( const nifti_image * const targetImageIn ) const;

  void CubicSplineTransformImageAdjoint3D(const nifti_image * const sourceImage, const nifti_image * const sourceWeightsImage, const nifti_image *warpedImage, const nifti_image * const warpedWeightsImage );
  void CubicSplineTransformImageAdjoint2D(const nifti_image * const sourceImage, const nifti_image * const sourceWeightsImage, const nifti_image *warpedImage, const nifti_image * const warpedWeightsImage );
  void NearestNeighbourTransformImageAdjoint3D( const nifti_image * const sourceImage, const nifti_image * const sourceWeightsImage, const nifti_image * warpedImage, const nifti_image * const warpedWeightsImage );
  void NearestNeighbourTransformImageAdjoint2D(const nifti_image * const sourceImage, const nifti_image *sourceWeightsImage, const nifti_image * const warpedImage, const nifti_image * const warpedWeightsImage );
  void TrilinearTransformImageAdjoint(const nifti_image * const sourceImage, const nifti_image * const sourceWeightsImage, const nifti_image *warpedImage, const nifti_image * const warpedWeightsImage );
  void BilinearTransformImageAdjoint(const nifti_image * const sourceImage, const nifti_image * const sourceWeightsImage, const nifti_image *warpedImage, const nifti_image * const warpedWeightsImage );

  bool dvfImageUpdateRequired;                 ///< Indicates if the DVF requires re-calculation when it is required
  unsigned int numberOfParameters;             ///< The total number of parameters needed to describe the transformation
  unsigned int numberOfLevelsToPerform;        ///< The total number of levels that are used to fit the model. Parameter needed to correctly scale constraints.
  int interpolation;                           ///< Type of interpolation used \todo Change into enumeration.
  nifti_image* deformationVectorFieldImage;    ///< The deformation vector field image
  PrecisionType warpedPaddingValue;            ///< Padding value for the warped image
  bool needToDeleteDVF;
};

}