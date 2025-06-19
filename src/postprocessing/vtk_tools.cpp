// #include "vtk_tools.hpp"


// vtk_tools::vtk_tools(){};


// void vtk_tools::write_vtk(const std::vector<std::vector<double>> &data,
//                    const std::string &vtk_name,
//                    double dx, double dt){
//     const int nt = static_cast<int>(data.size());
//     if (nt == 0) throw std::runtime_error("empty data");

//     const int nx = static_cast<int>(data[0].size());

//     // 1. 填充 VTK ImageData
//     auto img = vtkSmartPointer<vtkImageData>::New();
//     img->SetDimensions(nx, nt, 1);
//     img->SetOrigin(0.0, 0.0, 0.0);
//     img->SetSpacing(dx, dt, 1.0);

//     auto scal = vtkSmartPointer<vtkDoubleArray>::New();
//     scal->SetName("T");
//     scal->SetNumberOfComponents(1);
//     scal->SetNumberOfTuples(nx * nt);

//     // VTK ImageData 采用列主序 (i,j,k)
//     for (int j = 0; j < nt; ++j)
//         for (int i = 0; i < nx; ++i)
//             scal->SetTuple1(j * nx + i, data[j][i]);

//     img->GetPointData()->SetScalars(scal);

//     // 2. 写 .vti
//     auto wr = vtkSmartPointer<vtkXMLImageDataWriter>::New();
//     wr->SetFileName(vtkname.c_str());
//     wr->SetInputData(img);
//     wr->Write();
// }