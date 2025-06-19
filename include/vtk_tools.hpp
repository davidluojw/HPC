#ifndef VTK_TOOLS
#define VTK_TOOLS

#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <iostream>
#include <cmath>
#include <vector>

class vtk_tools
{
private:
    
public:
    vtk_tools();
    ~vtk_tools();

    // The function vtk_tools:: write_vtk is used to write two-dimensional spatiotemporal data 
    // (such as a one-dimensional temperature field that evolves over time) into a VTK format file (. vti) 
    // for visualization in tools such as ParaView. 
    // Its inputs include:
    // data: A two-dimensional vector that stores spatial data for multiple time steps (data [time steps] [spatial positions])
    // vtk_name: Output file name (such as "temperature. vti")
    // dx and dt: spatial step size (dx) and time step size (dt)
    void write_vtk(const std::vector<std::vector<double>> &data,
                   const std::string &vtk_name,
                   double dx, double dt, double zscale);


};




#endif