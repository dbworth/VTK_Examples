/*
Draw an annulus (disk or donut shape)
using 4 different methods:
 - vtkDiskSource
 - vtkContourTriangulator
 - vtkDelaunay2D
 - Lines

You can choose to draw the outline and/or the mesh,
and fill with a solid color or gradient.

Author:
David Butterworth, Oct 2018

Inspired by this sample code:
http://vtk.1045678.n5.nabble.com/How-to-set-the-boundary-of-vtkDelaunay2D-td5734142.html
*/

#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkColorTransferFunction.h>
#include <vtkContourTriangulator.h>
#include <vtkDelaunay2D.h>
#include <vtkDiskSource.h>
#include <vtkFeatureEdges.h>
#include <vtkLine.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkTextMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVector.h>
#include <vtkVertexGlyphFilter.h>

//#include <vtkButterflySubdivisionFilter.h>
#include <vtkLinearSubdivisionFilter.h>
//#include <vtkLoopSubdivisionFilter.h>


enum FillType
{
  None = 0,
  SolidColor,
  GradientColorRadius,
  GradientColorAngle
};

// Get the angle between two vectors in radians
const double angleBetweenVectors(const vtkVector3d& v1, const vtkVector3d& v2)
{
  const vtkVector3d v1n(v1.Normalized());
  const vtkVector3d v2n(v2.Normalized());

  const double a = v1n.Dot(v2n);
  double angle = std::acos(a);
  const vtkVector3d cross_product = v1n.Cross(v2n);

  // Get the angle from -M_PI to +M_PI
  const vtkVector3d normal(0.0, 0.0, 1.0);
  if (normal.Dot(cross_product) < 0.0)
  {
    angle = -1.0 * angle;
  }

  // Make all angles 0.0 to 2.0*M_PI
  if (angle < 0.0)
  {
    angle = 2.0*M_PI + angle;
  }

  return angle;
}

// Normalize an RGB value
void normalizeColor(double* rgb)
{
  // rgb[0] = Red
  // rgb[1] = Green
  // rgb[2] = Blue
  const double maxval = std::max(rgb[2], std::max(rgb[0], rgb[1]));
  rgb[0] = rgb[0] / maxval;
  rgb[1] = rgb[1] / maxval;
  rgb[2] = rgb[2] / maxval;
}

// Use a color transfer Function to generate the colors in the lookup table.
//
// ColorSpace = RGB by default
// Scale = linear by default
vtkSmartPointer<vtkLookupTable> makeLUTFromColorTransferFunction(vtkColorTransferFunction* ctf,
                                                                 const size_t lut_num_values,
                                                                 const double range_min, const double range_max)
{
  vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
  lut->SetNumberOfTableValues(lut_num_values);
  lut->SetTableRange(range_min, range_max);
  lut->Build();
 
  for (size_t i = 0; i < lut_num_values; ++i)
  {
    double* rgb;
    rgb = ctf->GetColor(static_cast<double>(i) / (lut_num_values-1));
    lut->SetTableValue(i, rgb);
  }

  return lut;
}

// Draw a disk with a hole removed
// using vtkDiskSource
void drawDisk(vtkActor* outlineActor, const bool drawOutline, double* outlineColor,
              vtkActor* meshActor, const bool drawMeshWireframe, double* meshColor,
              const int fillType, double* solidColor)
{
  // Parameters
  const double innerRadius = 1.0;
  const double outerRadius = 2.0;
  const int radialResolution = 2; // Number of points along the radius (default = 1)
  const int circumferentialResolution = 18; // Number of points around circumferance

  // Constants
  const double cx = 0.0; // Center of the disk. Don't need to change this
  const double cy = 0.0; // if using vtkActor->SetPosition(x,y,0)

  vtkSmartPointer<vtkDiskSource> diskSource = vtkSmartPointer<vtkDiskSource>::New();
  diskSource->SetInnerRadius(innerRadius); // default = 0.25
  diskSource->SetOuterRadius(outerRadius); // default = 0.5
  diskSource->SetRadialResolution(radialResolution); // default = 1
  diskSource->SetCircumferentialResolution(circumferentialResolution); // default = 6
  diskSource->Update();

  vtkSmartPointer<vtkPolyData> meshPolyData = vtkSmartPointer<vtkPolyData>::New();
  meshPolyData = diskSource->GetOutput();

  if (drawOutline)
  {
    vtkSmartPointer<vtkPolyDataMapper> outlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();

    // Extract edges from polygonal data (boundary, non-manifold, feature, manifold)
    vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
    #if VTK_MAJOR_VERSION <= 5
    featureEdges->SetInputConnection(diskSource->GetProducerPort());
    #else
    featureEdges->SetInputData(diskSource->GetOutput());
    #endif
    featureEdges->BoundaryEdgesOn();
    featureEdges->ManifoldEdgesOff();
    featureEdges->FeatureEdgesOff();
    featureEdges->NonManifoldEdgesOff();
    featureEdges->ColoringOff(); // Disable default colors e.g. boundary edges are red
    featureEdges->Update();

    #if VTK_MAJOR_VERSION <= 5
    outlineMapper->SetInputConnection(featureEdges->GetOutputPort());
    #else
    outlineMapper->SetInputData(featureEdges->GetOutput());
    #endif

    outlineActor->SetMapper(outlineMapper);
    outlineActor->GetProperty()->SetColor(outlineColor);
  }

  // If user doesn't want to render anything, then don't.
  // We need to check this because vtkDiskSource is
  // rendered as a solid white disk by default when
  // no other parameters are specified.
  if ((fillType == FillType::None) && (drawMeshWireframe == false))
  {
    return;
  }

  if (fillType == FillType::None)
  {
    if (drawMeshWireframe)
    {
      // Use solid color to draw edges of polyData
      meshActor->GetProperty()->SetRepresentationToWireframe();
      meshActor->GetProperty()->SetColor(meshColor);
    }
  }
  else if (fillType == FillType::SolidColor)
  {
    meshActor->GetProperty()->SetColor(solidColor);
  }
  else if ((fillType == FillType::GradientColorRadius) || (fillType == FillType::GradientColorAngle))
  {
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New(); // Defaults to RGB space
    vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New(); 

    // Set the known colors, the rest will be auto-generated.
    // Note: the function input is in the range 0.0 to 1.0
    if (fillType == FillType::GradientColorRadius)
    {
      // Linear gradient
      ctf->AddRGBPoint(0.0, 1.0, 0.0, 0.0); // Red
      ctf->AddRGBPoint(0.5, 0.0, 1.0, 0.0); // Green
      ctf->AddRGBPoint(1.0, 0.0, 0.0, 1.0); // Blue
      colorLookupTable = makeLUTFromColorTransferFunction(ctf, 3, innerRadius, outerRadius);
    }
    else if (fillType == FillType::GradientColorAngle)
    {
      // Circular gradient
      ctf->AddRGBPoint(0.0,         1.0, 0.0, 0.0); // Red
      ctf->AddRGBPoint(120.0/360.0, 0.0, 1.0, 0.0); // Green
      ctf->AddRGBPoint(240.0/360.0, 0.0, 0.0, 1.0); // Blue
      ctf->AddRGBPoint(1.0,         1.0, 0.0, 0.0); // Red
      colorLookupTable = makeLUTFromColorTransferFunction(ctf, 360, 0.0, 360.0);
    }

    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    // Generate a color for each point in the mesh based on the color map
    for (vtkIdType i = 0; i < meshPolyData->GetNumberOfPoints(); i++)
    {
      double p[3];
      meshPolyData->GetPoint(i, p);

      double val = 0.0;
      if (fillType == FillType::GradientColorRadius)
      {
        // Euclidean distance from mesh point to center of circle
        val = sqrt((p[0] - cx)*(p[0] - cx) + (p[1] - cy)*(p[1] - cy));          
      }

      else if (fillType == FillType::GradientColorAngle)
      {
        // Angle of rotation from 0 degrees (pointing up)
        const vtkVector3d v1(0.0,  1.0,  0.0);
        const vtkVector3d v2(p[0], p[1], 0.0);
        // Compute the angle such that the color gradient
        // is applied in a clockwise direction:
        val = vtkMath::DegreesFromRadians(angleBetweenVectors(v2, v1));
      }

      double dcolor[3];
      colorLookupTable->GetColor(val, dcolor); // Color components are scaled 0.0 to 1.0

      normalizeColor(dcolor);

      unsigned char color[3];
      for (unsigned int j = 0; j < 3; j++)
      {
        color[j] = std::round(255.0 * dcolor[j] / 1.0);
      }
      colors->InsertNextTupleValue(color); // Expects 0 to 255
    }

    meshPolyData->GetPointData()->SetScalars(colors);
  }

  vtkSmartPointer<vtkPolyDataMapper> meshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  #if VTK_MAJOR_VERSION <= 5
  meshMapper->SetInputConnection(meshPolyData->GetProducerPort());
  #else
  meshMapper->SetInputData(meshPolyData);
  #endif

  if ((fillType != FillType::None) && (drawMeshWireframe))
  {
    // Draw the wireframe on top of the fill color
    meshActor->GetProperty()->SetEdgeColor(meshColor);
    meshActor->GetProperty()->EdgeVisibilityOn();
    //meshActor->GetProperty()->SetInterpolationToFlat();
  }

  meshActor->SetMapper(meshMapper);
}

// Draw an annulus (2D donut shape)
// by generating two circles and meshing them
void drawAnnulus(vtkActor* outlineActor, const bool drawOutline, double* outlineColor,
                 vtkActor* meshActor, const bool drawMeshWireframe, double* meshColor,
                 const bool useDelaunayTriangulation, const bool enableMeshSubdivision,
                 const int fillType, double* solidColor)
{
  // Parameters
  const double innerRadius = 1.0;
  const double outerRadius = 2.0;
  const int circumferentialResolution = 18; // The annulus is a polygon of n sides
  const int mesh_number_of_subdivisions = 2;

  // Constants
  const double stepSize = 360.0 / static_cast<double>(circumferentialResolution);
  const double cx = 0.0; // Center of the annulus. Don't need to change this
  const double cy = 0.0; // if using vtkActor->SetPosition(x,y,0)
  const int num_steps = std::floor(360.0 / stepSize);

  vtkSmartPointer<vtkPoints> outlinePoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> outlineCellArray = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPolyData> outlinePolyData = vtkSmartPointer<vtkPolyData>::New();

  if (useDelaunayTriangulation)
  {
    // vtkDelaunay2D:

    outlineCellArray->InsertNextCell(num_steps);
    double theta = 0.0;
    while (theta  < 360.0)
    {
      const double x = outerRadius * cos(theta * M_PI/180.0);
      const double y = outerRadius * sin(theta * M_PI/180.0);
      vtkIdType pointId = outlinePoints->InsertNextPoint(x, y, 0.0);
      outlineCellArray->InsertCellPoint(pointId);
      theta += stepSize;
    }

    // Generate the inner circle (with negative winding order)
    // using line segments.
    outlineCellArray->InsertNextCell(num_steps);
    theta = 360.0 - stepSize;
    while (theta > -0.001) // theta >= 0.0
    {
      const double x = innerRadius * cos(theta * M_PI/180.0);
      const double y = innerRadius * sin(theta * M_PI/180.0);
      vtkIdType pointId = outlinePoints->InsertNextPoint(x, y, 0.0);
      outlineCellArray->InsertCellPoint(pointId);
      theta -= stepSize;
    }

    outlinePolyData->SetPoints(outlinePoints);
    outlinePolyData->SetPolys(outlineCellArray);
  }
  else
  {
    // vtkContourTriangulator:

    // Generate the outer circle (with positive winding order)
    // using line segments.
    vtkIdType prevPointId = num_steps - 1;
    double theta = 0.0;
    while (theta  < 360.0)
    {
      const double x = outerRadius * cos(theta * M_PI/180.0);
      const double y = outerRadius * sin(theta * M_PI/180.0);
      vtkIdType pointId = outlinePoints->InsertNextPoint(x, y, 0.0);
      outlineCellArray->InsertNextCell(2);
      outlineCellArray->InsertCellPoint(prevPointId);
      outlineCellArray->InsertCellPoint(pointId);
      prevPointId = pointId;
      theta += stepSize;
    }

    // Generate the inner circle (with negative winding order)
    // using line segments.
    prevPointId = (2 * num_steps) - 1;
    theta = 360.0 - stepSize;
    while (theta > -0.001) // theta >= 0.0
    {
      const double x = innerRadius * cos(theta * M_PI/180.0);
      const double y = innerRadius * sin(theta * M_PI/180.0);
      vtkIdType pointId = outlinePoints->InsertNextPoint(x, y, 0.0);
      outlineCellArray->InsertNextCell(2);
      outlineCellArray->InsertCellPoint(prevPointId);
      outlineCellArray->InsertCellPoint(pointId);
      prevPointId = pointId;
      theta -= stepSize;
    }

    outlinePolyData->SetPoints(outlinePoints);
    outlinePolyData->SetLines(outlineCellArray);
    outlinePolyData->BuildCells();
    outlinePolyData->BuildLinks();
  }

  if (drawOutline)
  {
    vtkSmartPointer<vtkPolyDataMapper> outlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();

    if (useDelaunayTriangulation)
    {
      // vtkDelaunay2D:

      // Extract edges from polygonal data (boundary, non-manifold, feature, manifold)
      vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
      #if VTK_MAJOR_VERSION <= 5
      featureEdges->SetInputConnection(outlinePolyData->GetProducerPort());
      #else
      featureEdges->SetInputData(outlinePolyData);
      #endif
      featureEdges->BoundaryEdgesOn();
      featureEdges->ManifoldEdgesOff();
      featureEdges->FeatureEdgesOff();
      featureEdges->NonManifoldEdgesOff();
      featureEdges->ColoringOff(); // Disable default colors e.g. boundary edges are red
      featureEdges->Update();

      #if VTK_MAJOR_VERSION <= 5
      outlineMapper->SetInputConnection(featureEdges->GetOutputPort());
      #else
      outlineMapper->SetInputData(featureEdges->GetOutput());
      #endif
    }
    else
    {
      // vtkContourTriangulator:

      #if VTK_MAJOR_VERSION <= 5
      outlineMapper->SetInputConnection(outlinePolyData->GetProducerPort());
      #else
      outlineMapper->SetInputData(outlinePolyData);
      #endif
    }

    outlineActor->SetMapper(outlineMapper);
    outlineActor->GetProperty()->SetColor(outlineColor);
  }

  if (drawMeshWireframe || (fillType != FillType::None))
  {
    //vtkSmartPointer<vtkPolyDataMapper> meshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkPolyData> meshPolyData = vtkSmartPointer<vtkPolyData>::New();
  
    if (useDelaunayTriangulation)
    {
      // vtkDelaunay2D:

      vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
      #if VTK_MAJOR_VERSION <= 5
      delaunay->SetInput(outlinePolyData);
      delaunay->SetSource(outlinePolyData);
      #else
      delaunay->SetInputData(outlinePolyData);
      delaunay->SetSourceData(outlinePolyData);
      #endif
      delaunay->Update();

      //meshMapper->SetInputConnection(delaunay->GetOutputPort());

      // Copy the PolyData, in case the user wants to add scalars
      meshPolyData = delaunay->GetOutput();
    }
    else
    {
      // vtkContourTriangulator:

      // Triangulate the points
      vtkSmartPointer<vtkContourTriangulator> triangulator = vtkSmartPointer<vtkContourTriangulator>::New();
      #if VTK_MAJOR_VERSION <= 5
      triangulator->SetInput(outlinePolyData);
      #else
      triangulator->SetInputData(outlinePolyData);
      #endif
      triangulator->Update();

      //meshMapper->SetInputConnection(triangulator->GetOutputPort());

      // Copy the PolyData, in case the user wants to add scalars
      meshPolyData = triangulator->GetOutput();
    }

    if (enableMeshSubdivision)
    {
      // Increase resolution of the mesh by sub-dividing each triangle.
      // For example, if n=2, the number of triangles in the resulting mesh
      // will be 16x the number of triangles in the original mesh.
      // This is required to set a smooth color gradient.
      vtkSmartPointer<vtkLinearSubdivisionFilter> subdivisionFilter = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
      subdivisionFilter->SetNumberOfSubdivisions(mesh_number_of_subdivisions);

      #if VTK_MAJOR_VERSION <= 5
      subdivisionFilter->SetInputConnection(meshPolyData->GetProducerPort());
      #else
      subdivisionFilter->SetInputData(meshPolyData);
      #endif
      subdivisionFilter->Update();

      meshPolyData = subdivisionFilter->GetOutput();
    }

    if (fillType == FillType::None)
    {
      // Use solid color to draw edges of polyData
      meshActor->GetProperty()->SetRepresentationToWireframe();
      meshActor->GetProperty()->SetColor(meshColor);
    }
    else if (fillType == FillType::SolidColor)
    {
      meshActor->GetProperty()->SetColor(solidColor);
    }
    else if ((fillType == FillType::GradientColorRadius) || (fillType == FillType::GradientColorAngle))
    {
      vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New(); // Defaults to RGB space
      vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New(); 

      // Set the known colors, the rest will be auto-generated.
      // Note: the function input is in the range 0.0 to 1.0
      if (fillType == FillType::GradientColorRadius)
      {
        // Linear gradient
        ctf->AddRGBPoint(0.0, 1.0, 0.0, 0.0); // Red
        ctf->AddRGBPoint(0.5, 0.0, 1.0, 0.0); // Green
        ctf->AddRGBPoint(1.0, 0.0, 0.0, 1.0); // Blue
        colorLookupTable = makeLUTFromColorTransferFunction(ctf, 3, innerRadius, outerRadius);
      }
      else if (fillType == FillType::GradientColorAngle)
      {
        // Circular gradient
        ctf->AddRGBPoint(0.0,         1.0, 0.0, 0.0); // Red
        ctf->AddRGBPoint(120.0/360.0, 0.0, 1.0, 0.0); // Green
        ctf->AddRGBPoint(240.0/360.0, 0.0, 0.0, 1.0); // Blue
        ctf->AddRGBPoint(1.0,         1.0, 0.0, 0.0); // Red
        colorLookupTable = makeLUTFromColorTransferFunction(ctf, 360, 0.0, 360.0);
      }

      vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
      colors->SetNumberOfComponents(3);
      colors->SetName("Colors");

      // Generate a color for each point in the mesh based on the color map
      for(vtkIdType i = 0; i < meshPolyData->GetNumberOfPoints(); i++)
      {
        double p[3];
        meshPolyData->GetPoint(i, p);

        double val = 0.0;
        if (fillType == FillType::GradientColorRadius)
        {
          // Euclidean distance from mesh point to center of circle
          val = sqrt((p[0] - cx)*(p[0] - cx) + (p[1] - cy)*(p[1] - cy));          
        }

        else if (fillType == FillType::GradientColorAngle)
        {
          // Angle of rotation from 0 degrees (pointing up)
          const vtkVector3d v1(0.0,  1.0,  0.0);
          const vtkVector3d v2(p[0], p[1], 0.0);
          // Compute the angle such that the color gradient
          // is applied in a clockwise direction:
          val = vtkMath::DegreesFromRadians(angleBetweenVectors(v2, v1));
        }

        double dcolor[3];
        colorLookupTable->GetColor(val, dcolor); // Color components are scaled 0.0 to 1.0

        normalizeColor(dcolor);

        unsigned char color[3];
        for (unsigned int j = 0; j < 3; j++)
        {
          color[j] = std::round(255.0 * dcolor[j] / 1.0);
        }
        colors->InsertNextTupleValue(color); // Expects 0 to 255
      }

      meshPolyData->GetPointData()->SetScalars(colors);
    }

    vtkSmartPointer<vtkPolyDataMapper> meshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    #if VTK_MAJOR_VERSION <= 5
    meshMapper->SetInputConnection(meshPolyData->GetProducerPort());
    #else
    meshMapper->SetInputData(meshPolyData);
    #endif

    if ((fillType != FillType::None) && (drawMeshWireframe))
    {
      // Draw the wireframe on top of the fill color
      meshActor->GetProperty()->SetEdgeColor(meshColor);
      meshActor->GetProperty()->EdgeVisibilityOn();
      //meshActor->GetProperty()->SetInterpolationToFlat();
    }

    meshActor->SetMapper(meshMapper);
  }
}

// Draw a wheel shape
// using lines along the radius
void drawWheel(vtkActor* endpointsActor, const bool drawEndPoints, double* endpointsColor,
               vtkActor* linesActor,
               const int fillType, double* solidColor)
{
  // Parameters
  const double innerRadius = 1.0;
  const double outerRadius = 2.0;
  const int circumferentialResolution = 800; // Draw a wheel using n lines e.g. 800 ==> 0.0025 (radians)
  const int lineWidth = 2; // default = 1

  // Constants
  const double cx = 0.0; // Center of the wheel. Don't need to change this
  const double cy = 0.0; // if using vtkActor->SetPosition(x,y,0)
  const double step_size = 2.0 / static_cast<double>(circumferentialResolution); // num steps from -1.0 to +1.0

  vtkSmartPointer<vtkPoints> linePoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> linesCellArray = vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);

  vtkIdType ptIdx = 0;

  // Draw a wheel by drawing individual lines
  // along the radius, from angle -PI to +PI:

  // Iterate from -1.0 to 0.99
  // i = arctan(Y,X)/PI
  unsigned int count = 0;
  double i = -1.0;
  while (i < 1.0)
  {
    // Compute the line coordinates
    const double x1 = innerRadius * cos(i * M_PI);
    const double y1 = innerRadius * sin(i * M_PI);
    const double x2 = outerRadius * cos(i * M_PI);
    const double y2 = outerRadius * sin(i * M_PI);

    i = i + step_size;

    // Double-check if we have drawn all the lines.
    // Sometimes the while() condition fails!
    // e.g.
    // 1.0667 < 1.0
    // 1.0025 < 1.0
    count++;
    if (count > circumferentialResolution)
    {
      break;
    }

    // Add two points
    const double p0[3] = {cx + x1,
                          cy + y1,
                          0.0};
    const double p1[3] = {cx + x2,
                          cy + y2,
                          0.0};
    linePoints->InsertNextPoint(p0);
    linePoints->InsertNextPoint(p1);

    if (fillType != FillType::None)
    {
      // Add a line
      vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
      line->GetPointIds()->SetId(0, ptIdx); // (0, Index of 'p0' in 'linePoints')
      line->GetPointIds()->SetId(1, ptIdx+1); // (1, Index of 'p1' in 'linePoints')
      linesCellArray->InsertNextCell(line);
      ptIdx += 2;
    }

    if (fillType == FillType::GradientColorAngle)
    {
      // Generate the (R,G,B) colors:
      //   R = 1 at A = 0
      //   G = 1 at A = 2/3
      //   B = 1 at A = -2/3 

      double R = 1.0 - std::fabs(i); // R goes from 0.0 to 1.0 to 0.01

      double AA = i - 2.0/3.0;
      if (AA < -1.0)
      {
        AA = 2.0 + AA;
      }
      double G = 1.0 - std::fabs(AA);

      AA = i + 2.0/3.0;
      if (AA > 1.0)
      {
        AA = 2.0 - AA;
      }
      double B = 1.0 - std::fabs(AA);

      // Boost low end (brighten colors)   
      //R = R + 0.2
      //G = G + 0.2
      //B = B + 0.2

      double dcolor[3] = {R, G, B};

      normalizeColor(dcolor);

      unsigned char color[3];
      for (unsigned int j = 0; j < 3; j++)
      {
        color[j] = std::round(255.0 * dcolor[j] / 1.0);
      }
      colors->InsertNextTupleValue(color); // Expects 0 to 255
    }
  }

  vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
  linesPolyData->SetPoints(linePoints);

  if (drawEndPoints)
  {
    // To draw the end-points of each line, we need to
    // filter the points through vtkVertexGlyphFilter:

    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
    vertexFilter->SetInputConnection(linesPolyData->GetProducerPort());
    #else
    vertexFilter->SetInputData(linesPolyData);
    #endif
    vertexFilter->Update();

    vtkSmartPointer<vtkPolyData> endpointsPolyData = vtkSmartPointer<vtkPolyData>::New();
    endpointsPolyData->ShallowCopy(vertexFilter->GetOutput());
   
    vtkSmartPointer<vtkPolyDataMapper> endpointsMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
   
    #if VTK_MAJOR_VERSION <= 5
    endpointsMapper->SetInput(endpointsPolyData);
    #else
    endpointsMapper->SetInputData(endpointsPolyData);
    #endif

    endpointsActor->SetMapper(endpointsMapper);
    endpointsActor->GetProperty()->SetPointSize(1);
    endpointsActor->GetProperty()->SetColor(endpointsColor);
  }

  if (fillType == FillType::None)
  {
    // Not supported, don't render anything
    return;
  }
  else if (fillType == FillType::SolidColor)
  {
    linesActor->GetProperty()->SetColor(solidColor);
  }
  else if (fillType == FillType::GradientColorRadius)
  {
    // Not supported, don't render anything
    return;
  }
  else if (fillType == FillType::GradientColorAngle)
  {
    // Color the lines
    linesPolyData->GetCellData()->SetScalars(colors);
  }

  linesPolyData->SetLines(linesCellArray);

  vtkSmartPointer<vtkPolyDataMapper> linesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();

  #if VTK_MAJOR_VERSION <= 5
  linesMapper->SetInput(linesPolyData);
  #else
  linesMapper->SetInputData(linesPolyData);
  #endif

  linesActor->SetMapper(linesMapper);
  linesActor->GetProperty()->SetLineWidth(lineWidth);
}


int main(int, char *[])
{
  // Options
  const bool drawMeshWireframe = true;
  double meshColor[3] = {0, 0, 1}; // blue

  const int fillType = FillType::SolidColor; // Options: None, SolidColor, GradientColorRadius, GradientColorAngle
  double solidColor[3] = {1.0, 0.650, 0.298}; // light orange

  // Options only for drawDisk(), drawAnnulus() and drawWheel()
  const bool drawOutline = true;
  //double outlineColor[3] = {1, 0, 0}; // red
  double outlineColor[3] = {1, 1, 1}; // white

  // Options only for drawAnnulus()
  const bool enableMeshSubdivision = true; // True = increase number of triangles by 16x


  // Render an annulus using vtkDiskSource (the most efficient)
  vtkSmartPointer<vtkActor> diskMeshActor = vtkSmartPointer<vtkActor>::New();
  vtkSmartPointer<vtkActor> diskOutlineActor = vtkSmartPointer<vtkActor>::New();
  drawDisk(diskOutlineActor, drawOutline, outlineColor,
           diskMeshActor, drawMeshWireframe, meshColor,
           fillType, solidColor);

  // Render an annulus by meshing two circles with Contour Triangulator
  bool useDelaunayTriangulation = false; // False = vtkContourTriangulator
  vtkSmartPointer<vtkActor> annulusOutlineActor1 = vtkSmartPointer<vtkActor>::New();
  vtkSmartPointer<vtkActor> annulusMeshActor1 = vtkSmartPointer<vtkActor>::New();
  drawAnnulus(annulusOutlineActor1, drawOutline, outlineColor,
              annulusMeshActor1, drawMeshWireframe, meshColor,
              useDelaunayTriangulation, enableMeshSubdivision,
              fillType, solidColor);

  // Render an annulus by meshing two circles with Delaunay triangulation
  useDelaunayTriangulation = true; // True = vtkDelaunay2D
  vtkSmartPointer<vtkActor> annulusOutlineActor2 = vtkSmartPointer<vtkActor>::New();
  vtkSmartPointer<vtkActor> annulusMeshActor2 = vtkSmartPointer<vtkActor>::New();
  drawAnnulus(annulusOutlineActor2, drawOutline, outlineColor,
              annulusMeshActor2, drawMeshWireframe, meshColor,
              useDelaunayTriangulation, enableMeshSubdivision,
              fillType, solidColor);

  // Render an annulus by drawing multiple lines along the radius
  vtkSmartPointer<vtkActor> wheelEndpointsActor = vtkSmartPointer<vtkActor>::New();
  vtkSmartPointer<vtkActor> wheelLinesActor = vtkSmartPointer<vtkActor>::New();
  drawWheel(wheelEndpointsActor, drawOutline, outlineColor,
            wheelLinesActor,
            fillType, solidColor);


  std::vector<std::string> names;
  names = std::vector<std::string>({"Mesh (vtkDiskSource)",
                                    "Mesh (vtkContourTriangulator)",
                                    "Mesh (vtkDelaunay2D)",
                                    "Lines"});

  std::vector<vtkSmartPointer<vtkActor> > actors = {diskMeshActor,
                                                    annulusMeshActor1,
                                                    annulusMeshActor2,
                                                    wheelLinesActor
                                                   };

  // Create the labels
  std::vector<vtkSmartPointer<vtkActor2D> > textActors;
  for (unsigned int i = 0; i < actors.size(); i++)
  {
    vtkSmartPointer<vtkTextMapper> textMapper = vtkSmartPointer<vtkTextMapper>::New();
    //textMapper->SetInput(actors[i]->GetClassName()); // set text label to the name of the object source
    textMapper->SetInput(names[i].c_str());
 
    textActors.push_back(vtkSmartPointer<vtkActor2D>::New());
    textActors[i]->SetMapper(textMapper);
    textActors[i]->SetPosition(10, 10); // Specified in display coordinates
  }

  std::vector<vtkSmartPointer<vtkActor> > outlineActors = {diskOutlineActor,
                                                           annulusOutlineActor1,
                                                           annulusOutlineActor2,
                                                           wheelEndpointsActor
                                                          };



  /*
  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // Add the actors to the scene
  if (drawMeshWireframe || (fillType != FillType::None)) // fill_solid)
  {
    //renderer->AddActor(annulusMeshActor);
  }
  //renderer->AddActor(annulusOutlineActor);

  //renderer->AddActor(wheelActor);

  renderer->AddActor(diskActor);

  renderer->SetBackground(48.0/255, 10.0/255, 36.0/255); // off-black

  // Render an image (lights and cameras are created automatically)
  renderWindow->SetWindowName("test");
  renderWindow->SetSize(640, 480); //(width, height)
  renderWindow->Render();

  // Top-down image view. Can only drag and zoom.
  vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
  renderWindowInteractor->SetInteractorStyle(style);

  renderWindowInteractor->Start(); // start the event loop
  */




  // Define size of the grid that will hold the objects
  int gridCols = 2;
  int gridRows = 2;

  // Define side length (in pixels) of each renderer square
  int rendererSize = 200;
 
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetSize(rendererSize*gridCols, rendererSize*gridRows);
 
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
 
  // Set up a grid of viewports for each renderer
  for (int row = 0; row < gridRows; row++)
  {
    for (int col = 0; col < gridCols; col++)
    {
      const double index = row*gridCols + col;
 
      // Create a renderer for this grid cell
      vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
      renderer->SetBackground(.6, .6, .6);
 
      // Set the renderer's viewport dimensions (xmin, ymin, xmax, ymax)
      // within the render window
      double viewport[4] = {
        static_cast<double>((col    )              * rendererSize) / (gridCols * rendererSize),
        static_cast<double>((gridRows - (row + 1)) * rendererSize) / (gridRows * rendererSize),
        static_cast<double>((col + 1)              * rendererSize) / (gridCols * rendererSize),
        static_cast<double>((gridRows - (row    )) * rendererSize) / (gridRows * rendererSize) };
      renderer->SetViewport(viewport);
 
      // Add the corresponding actor and label for this grid cell, if they exist
      if (index < actors.size())
      {
        renderer->AddActor(textActors[index]);
        renderer->AddActor(actors[index]);
        renderer->AddActor(outlineActors[index]);
      }




      /*
      if (index == 0)
      {
        renderer->AddActor(diskOutlineActor);
      }
      else if (index == 1)
      {
        renderer->AddActor(annulusOutlineActor1);
      }
      else if (index == 2)
      {
        renderer->AddActor(annulusOutlineActor2);
      }
      else if (index == 3)
      {
        renderer->AddActor(wheelEndpointsActor);
      }
      */

      renderWindow->AddRenderer(renderer);
    }
  }
 
  vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(renderWindow);
 
  renderWindow->Render();
  interactor->Start();

  return EXIT_SUCCESS;
}
