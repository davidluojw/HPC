#include "vtk_tools.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"

vtk_tools::vtk_tools(){}

vtk_tools::~vtk_tools() {}


void vtk_tools::write_vtk(const std::vector<std::vector<double>> &data,
                          const std::string &vtk_name,
                          double dx, double dt, double zscale)
{
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // only read h5 file on main process
    if (rank == 0){
        const int nt = static_cast<int>(data.size());         // Time steps
        std::cout << "nt: " << nt << std::endl;
        if (nt == 0) std::cerr << "Error creating dataspace" << std::endl;  // 

        const int nx = static_cast<int>(data[0].size());     // Space points 
        std::cout << "nx: " << nx << std::endl;

        const int np = nt * nx;  

        // 1. Fill VTK ImageData
        // vtkImageData: is used to represent regular grid data (such as spatiotemporal fields).
        // SetDimensions(nx, nt, 1): maps spatiotemporal data into a two-dimensional image (X-axis: spatial position, Y-axis: time step, Z-axis fixed at 1).
        // SetSpacing(dx, dt, 1.0): Set the grid spacing to be dx in the X direction, dt in the Y direction, and 1.0 in the Z direction (without practical significance)
        auto img = vtkSmartPointer<vtkImageData>::New();  
        img->SetDimensions(nx, nt, 1);       // Dimension: number of spatial points x number of time steps x 1 (two-dimensional)
        img->SetOrigin(0.0, 0.0, 0.0);       // Origin of coordinate system
        img->SetSpacing(dx, dt, 1.0);        // Grid spacing: dx (space), dt (time), 1 (meaningless)

        // vtkDoubleArray: stores double precision floating-point data,
        // SetName("T"):  named "T" (such as temperature).
        auto points = vtkSmartPointer<vtkPoints>::New();
        auto scal = vtkSmartPointer<vtkDoubleArray>::New();
        scal->SetName("T");                   // Field name (displayed as variable name in ParaView)
        scal->SetNumberOfComponents(1);       // Scalar data (single component)
        scal->SetNumberOfTuples(nx * nt);     // Total number of data points=number of spatial points x number of time steps
        std::cout << "nx: " << nx << std::endl;

        // VTK ImageData uses column main order (i, j, k)
        // Fill data in column main order (VTK storage order: traverse X first, then Y, and finally Z)
        // Column main sequence storage: The memory layout is (x0, t0) → (x1, t0) → .. → (x_{nx-1}, t0) → (x0, t1) → . ... complies with VTK rules.
        // Each grid point is associated with a scalar value (temperature), which is bound to vtkImageData through tScalars
        for (int j = 0; j < nt; ++j){                     // j: Time step index (Y-axis)
            for (int i = 0; i < nx; ++i){                 // i: Spatial Position Index (X-axis)
                const double x = i * dx;
                const double t = j * dt;
                const double z =  data[j][i] * zscale;    // 
                points->InsertNextPoint(x, t, z);
                scal->SetTuple1(j * nx + i, data[j][i]);   // Index=j * nx+i
            }
        }

        auto cells = vtkSmartPointer<vtkCellArray>::New();
        // 插入单元：每个单元是一个四边形，包含 4 个点
        for (int j = 0; j < nt - 1; ++j)
        {
            for (int i = 0; i < nx - 1; ++i)
            {
                vtkIdType id0 = j * nx + i;
                vtkIdType id1 = j * nx + i + 1;
                vtkIdType id2 = (j + 1) * nx + i + 1;
                vtkIdType id3 = (j + 1) * nx + i;

                auto quad = vtkSmartPointer<vtkQuad>::New();
                quad->GetPointIds()->SetId(0, id0);
                quad->GetPointIds()->SetId(1, id1);
                quad->GetPointIds()->SetId(2, id2);
                quad->GetPointIds()->SetId(3, id3);

                cells->InsertNextCell(quad);
            }
        }

        auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        grid->SetPoints(points);
        grid->SetCells(VTK_QUAD, cells);
        grid->GetPointData()->SetScalars(scal);

        img->GetPointData()->SetScalars(scal);             // Attach scalar data to grid points

        // 2. Write .vti
        // vtkXMLImageDataWriter: outputs VTkImageData as an XML formatted. vti file (VTK ImageData).
        // The process conforms to the VTK Writer standard steps: instantiation → setting file name → binding data → calling Write()
        
        auto wr = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        std::string filename = vtk_name + ".vti";
        wr->SetFileName(filename.c_str());   // Set output file name
        wr->SetInputData(img);               // Input VTK data object
        wr->Write();                         // Perform write operation

        auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        filename = vtk_name + ".vtu";
        writer->SetFileName(filename.c_str() );
        writer->SetInputData(grid);
        // writer->Write();
        try {
            writer->Write();
        } catch (const std::exception& e) {
            std::cerr << "VTK write error: " << e.what() << std::endl;
        }

        std::cout << "Done: " << vtk_name << "\n";
    }
}