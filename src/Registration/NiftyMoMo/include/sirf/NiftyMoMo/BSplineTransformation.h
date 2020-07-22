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

#include "_reg_localTrans.h"
#include "_reg_resampling.h"
#include "nifti1_io.h"
#include "Transformation.h"

namespace NiftyMoMo {

template<class DTYPE>
inline void get_BSplineBasisValues( DTYPE basis, DTYPE *values )
{
  DTYPE FF = basis*basis;
  DTYPE FFF = FF*basis;
  DTYPE MF = static_cast<DTYPE>(1.0 - basis);
  values[0] = static_cast<DTYPE>((MF) *(MF) *(MF) / (6.0));
  values[1] = static_cast<DTYPE>((3.0*FFF - 6.0*FF + 4.0) / 6.0);
  values[2] = static_cast<DTYPE>((-3.0*FFF + 3.0*FF + 3.0*basis + 1.0) / 6.0);
  values[3] = static_cast<DTYPE>(FFF / 6.0);
}


/** Class implementing the standard b-spline transformation 
 */
class BSplineTransformation :public Transformation
{
public:
  typedef Transformation::PrecisionType PrecisionType;

  /** Constructor.
   *  Requires the reference image to define the space in which the transformation is valid, the number of levels 
   *  to update the bending energy correctly as well as the final control-point grid spacing.
   *  \param referenceImageIn Pointer to the reference image
   *  \param numberOfLevelsToPerform Number of levels to perform
   *  \param finalControlPointGridSpacing Pointer to the final control-point grid spacing if negative, the spacing is assumed to be given 
   *                                      as multiples of the voxel spacing
   */
  BSplineTransformation( nifti_image* referenceImageIn, 
                         unsigned int numberOfLevelsToPerform, 
                         float* finalControlPointGridSpacing );
  
  /** Copy constructor.
   *  \param transformToCopy Reference to the B-spline transformation object to copy
   */
  BSplineTransformation( const BSplineTransformation& transformToCopy );


  /** Destructor.
   */
  ~BSplineTransformation();
  
  /** Initialise a given level, i.e. perform a grid refinement.
   *  If the same level was initialised before, no refinement will be performed. Currently it is not checked, if 
   *  a coarser level is requested.
   *  \param level Defines the level that is to be initialised.
   */
  virtual void InitialiseLevel( unsigned int level );
  
  /** Set the parameters that define the b-spline transformation.
   *  \param paramsIn Pointer to the parameters that will be copied to the internal transformation parameters.
   *  \param parametersAreDisplacements Set to true if the parameters are displacements (i.e. relative to control-point position and NOT absolute positions)
   */
  void SetParameters( PrecisionType* paramsIn, bool parametersAreDisplacements );
  
  /** For a given target image, calculate a deformation field. 
   *  \param targetImageIn The target image defines the geometry of the DVF, or the voxel locations for which deformation vectors are calculated.
   */
  virtual nifti_image* GetDeformationVectorField(const nifti_image * const targetImageIn );
   
  /** Calculates the gradient of an image for a given deformation vector field.
   *  \param denseDVFIn Deformation vector field 
   *  \param sourceImage The image for which the gradient is calculated
   */
  virtual BSplineTransformation::PrecisionType* GetDVFGradientWRTTransformationParameters( nifti_image* denseDVFIn );
  
  /** Calculate the gradient of the constraint term (regularisation) for the transformation with the current parameters.
   */
  virtual BSplineTransformation::PrecisionType* GetConstraintGradientWRTTransformationParameters();
  
  /** Get the constraint value for the transformation with the current parametes.
   */
  double GetConstraintValue();
  
  /** Return a pointer to the parameters of the transformation
   */
  virtual PrecisionType* GetParameters();
  
  /** Get the total number of transformation parameters (for the current level)
   */
  virtual unsigned int GetNumberOfParameters();
  
  /** B-spline specific regularisation term weight. Penalises the linear energy.
   *  \param linearEnergyWeightIn Input weight
   */
  void SetLinearEnergyWeight( double linearEnergyWeightIn );
  
  /** B-spline specific regularisation term weight. Penalises the bending energy.
   *  \param bendingEnergyWeightIn Input weight
   */
  void SetBendingEnergyWeight( double bendingEnergyWeightIn );
  
  /** Function to create a deep copy of this object.
   */
  std::shared_ptr<Transformation> DeepCopy();
  
  /** Display function to print the transformation parameters to the standard output
   */
  void DisplayTransformationParameters();
  
  /** Calculates the sum of the penalty term weights. Curretnly the bending energy 
   *  weight and the linear energy weight are added. Responsibility of the correct 
   *  scaling of the penalty term lies with the generating function.
   */
  virtual PrecisionType GetSumOfPenaltyWeights();

  /** Calculate the maximum control-point displacement value for the given transformation parametes. 
   */
  virtual PrecisionType GetMaxTransformationParameterLength( PrecisionType* parametersIn );


  /** Returns the b-spline control-point grind as an image. 
   */
  nifti_image* GetTransformationAsImage();

protected:

  /** Compute a dense vector field (DVF) from a B-spline control point grid (CPG).
   *  This function assumes the axes of the DVF is aligned with axes of the CPG
   *  \param denseVectorFieldImage: A nifti image that will be filled with the DVF
   */
  void CPGToDVF2D( nifti_image *denseVectorFieldImage );

  /** Compute a dense vector field (DVF) from a B-spline control point grid (CPG).
   *  This function assumes the axes of the DVF is aligned with axes of the CPG
   *  \param denseVectorFieldImage: A nifti image that will be filled with the DVF
   */
  void CPGToDVF3D( nifti_image *denseVectorFieldImage );

  /** Compute a B-spline control point grid (CPG) from the B-spline convolution
   *  of a dense 2D vector field (DVF), i.e. the CPG will be the B-spline parameterisation that best
   *  approximates the DVF. This function assumes the axes of the DVF are aligned with axes of the CPG
   *  \param controlPointGridImage: A nifti image that will be filled with the B-spline CPG
   *  \param denseVectorFieldImage: A nifti image containing the DVF
   *  \param normaliseKernel: If true each B-spline kernel will be normalised to have a sum of 1
   */
  void DVFToCPG2D( nifti_image *controlPointGridImage,
                   nifti_image *denseVectorFieldImage,
                   bool normaliseKernel );
  
  /** Compute a B-spline control point grid (CPG) from the B-spline convolution
   *  of a dense 3D vector field (DVF), i.e. the CPG will be the B-spline parameterisation that best
   *  approximates the DVF. This function assumes the axes of the DVF are aligned with axes of the CPG
   *  \param controlPointGridImage: A nifti image that will be filled with the B-spline CPG
   *  \param denseVectorFieldImage: A nifti image containing the DVF
   *  \param normaliseKernel: If true each B-spline kernel will be normalised to have a sum of 1
   */
  void DVFToCPG3D( nifti_image *controlPointGridImage,
                   nifti_image *denseVectorFieldImage,
                   bool normaliseKernel );


  ///** Function to reorientate vectors in 5D nifti image. 
  // *  Vectors which are stored along the 5th (u) dimension will be reoriented
  // *  \param vectorFieldToReorientate Vector image to be reorientated
  // *  \param reorientationMatrix 4 x 4 Matrix used to reorientate vectors
  // */
  //void ReorientateVectorImage( nifti_image* vectorFieldToReorientate, mat44 reorientationMatrix );

  nifti_image* controlPointGridImage;          ///< The control-point grid image
  double bendingEnergyWeight;                  ///< The bending-energy weight
  double linearEnergyWeight;                   ///< The linear energy weight

private:
  int lastInitialisedLevel;                    ///< The level taht was last requested to be initialised
};

}
