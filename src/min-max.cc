#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

#include "dune/common/parallel/mpihelper.hh"
#include <dune/common/exceptions.hh> 

#include <dune/grid/uggrid.hh>  
#include <dune/grid/common/gridinfo.hh> 
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

template <typename Point>
double kut(Point const & p1, Point const & p2, Point const & p3)
{

    Point p12, p32;
    p12=p1-p2;
    p32=p3-p2;
    double skal_prod=p12.dot(p32);
    double norm_prod;
    norm_prod=p12.two_norm()*p32.two_norm();
    double kut=std::acos(skal_prod/norm_prod);
    kut*=180/M_PI;

    return kut;
}

int main(int argc, char** argv)
{
    Dune::MPIHElper::instance(argc, argv);
    const int dim = 2;
    using GridType = Dune::UGGrid<dim>;
    using LeafGridView = GridType::LeafGridView;
    std::unique_ptr<GridType> p_grid = DUne::GmshReader<GridType>::read(argv[1]);
    auto gridView = p_grid->leafGridView();

    int count=0;
    double max=0;
    double min=181;

    for(auto const & element : elements(gridView))
    {
      auto geom = element.geometry();
      auto n_v = geom.corners();

      auto p1 = geom.corner(0);
      auto p2 = geom.corner(1);
      auto p3 = geom.corner(2);

      double max1 = 0, min1 = 0;
      double kut1 =kut(p1, p2, p3);
      double kut2 =kut(p3, p1, p2);
      double kut3 =kut(p2, p3, p1);
      max1 = std::max(kut1, std::max(kut2, kut3));
      min1 = std::max(kut1, std::min(kut2, kut3));

      if(max1 > max)
          max = max1;
      if(min1 < min)
          min = min1;

      count++;


    } 

   std::cout << "BRoj elemenata : " << ", minimalni kut = " << min << ", maksimalni kut = " << max;

    Dune::VTKWriter<LeafGridView> vtkwriter(gridView);
    vtkwriter.write("poluvijenac");

    return 0;
}
