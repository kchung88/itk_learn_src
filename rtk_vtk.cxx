// RTK includes
#include "rtkConfiguration.h"
#include <rtkFDKBackProjectionImageFilter.h>
#include <rtkConstantImageSource.h>
#include <rtkThreeDCircularProjectionGeometry.h>
#include <rtkRayEllipsoidIntersectionImageFilter.h>
#include <rtkDisplacedDetectorImageFilter.h>
#include <rtkParkerShortScanImageFilter.h>
#include <rtkFDKConeBeamReconstructionFilter.h>
#include <rtkCudaFDKConeBeamReconstructionFilter.h>

// ITK includes
#include <itkImageFileWriter.h>
#include <itkStreamingImageFilter.h>
#include "itkCudaImage.h"
#include "itkViewImage.h"




int main(int argc, char **argv)
{
	const unsigned int Dimension = 3;
	typedef float      OutputPixelType;
	// Defines the image type
	#ifdef RTK_USE_CUDA
		typedef itk::CudaImage< OutputPixelType, Dimension > ImageType;
	#else
		typedef itk::Image< OutputPixelType, Dimension > ImageType;
	#endif

	// Defines the RTK geometry object
	typedef rtk::ThreeDCircularProjectionGeometry GeometryType;
	GeometryType::Pointer geometry = GeometryType::New();

	// Projection matrices
	unsigned int numberOfProjections = 150;
	unsigned int firstAngle = 0;
	unsigned int angularArc = 360;
	unsigned int sid = 600; // source to isocenter distance in mm
	unsigned int sdd = 1200; // source to detector distance in mm
	int isox = 0; // X coordinate on the projection image of isocenter
	int isoy = 0; // Y coordinate on the projection image of isocenter

	for (unsigned int noProj = 0; noProj < numberOfProjections; noProj++)
	{
		double angle = (float)firstAngle + (float)noProj * angularArc / (float)numberOfProjections;
		geometry->AddProjection(sid,
			sdd,
			angle,
			isox,
			isoy);
	}

	// Create a stack of empty projection images
	typedef rtk::ConstantImageSource< ImageType > ConstantImageSourceType;
	ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
	ConstantImageSourceType::PointType origin;
	ConstantImageSourceType::SpacingType spacing;
	ConstantImageSourceType::SizeType sizeOutput;

	origin[0] = -129.;
	origin[1] = -129.;
	origin[2] = -129.;

	// Adjust size according to geometry
	sizeOutput[0] = 256;
	sizeOutput[1] = 256;
	sizeOutput[2] = numberOfProjections;

	spacing[0] = 1.;
	spacing[1] = 1.;
	spacing[2] = 1.;
	constantImageSource->SetOrigin(origin);
	constantImageSource->SetSpacing(spacing);
	constantImageSource->SetSize(sizeOutput);
	constantImageSource->SetConstant(0.);

	std::cout << "Setup done\n";
	std::cout << "Deine Mutter \n";

	// Create the projector
	typedef rtk::RayEllipsoidIntersectionImageFilter<ImageType, ImageType> REIType;
	REIType::Pointer rei = REIType::New();
	REIType::VectorType semiprincipalaxis, center;
	semiprincipalaxis.Fill(50.);
	center.Fill(0.);
	//Set GrayScale value, axes, center...
	rei->SetDensity(2.);
	rei->SetAngle(0.);
	rei->SetCenter(center);
	rei->SetAxis(semiprincipalaxis);
	rei->SetGeometry(geometry);
	rei->SetInput(constantImageSource->GetOutput());
	rei->Update();

	std::cout << "DRR done\n";

	using viewimage = itk::ViewImage< ImageType >;
	viewimage::View(rei->GetOutput());

	// Create reconstructed image
	ConstantImageSourceType::Pointer constantImageSource2 = ConstantImageSourceType::New();

	// Adjust size according to geometry
	sizeOutput[0] = 256;
	sizeOutput[1] = 256;
	sizeOutput[2] = 256;

	constantImageSource2->SetOrigin(origin);
	constantImageSource2->SetSpacing(spacing);
	constantImageSource2->SetSize(sizeOutput);
	constantImageSource2->SetConstant(0.);

	// FDK reconstruction filtering
	typedef rtk::CudaFDKConeBeamReconstructionFilter FDKGPUType;
	//typedef rtk::FDKConeBeamReconstructionFilter< ImageType > FDKCPUType;
	FDKGPUType::Pointer feldkamp = FDKGPUType::New();
	feldkamp->SetInput(0, constantImageSource2->GetOutput());
	feldkamp->SetInput(1, rei->GetOutput());
	feldkamp->SetGeometry(geometry);
	feldkamp->GetRampFilter()->SetTruncationCorrection(0.);
	feldkamp->GetRampFilter()->SetHannCutFrequency(0.0);
	feldkamp->Update();

	std::cout << "Reconstruction done\n";

	// Writer
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("output.tiff");
	writer->SetUseCompression(true);
	writer->SetInput(feldkamp->GetOutput());
	writer->Update();

	writer->SetFileName("proj.tiff");
	writer->SetUseCompression(true);
	writer->SetInput(rei->GetOutput());
	writer->Update();

	std::cout << "Writing done\n";
	return 0;
}