//
//=======================================================================
// Copyright 2016
// Author: Alex Hagen-Zanker
// University of Surrey
//
// Distributed under the MIT Licence (http://opensource.org/licenses/MIT)
//=======================================================================

#include <blink/raster_tools/distance_transform.h>
#include <blink/raster_tools/fuzzy_kappa.h>
#include <blink/raster/utility.h>

#include <blink/iterator/zip_range.h>

#include <ctime>

int demo_distance()
{
  auto input = blink::raster::open_gdal_raster<int>("input.tif", GA_ReadOnly);
  auto output_euclidean = blink::raster::create_gdal_raster_from_model<double>(
    "output_euclidean.tif", input);
  auto output_squared_euclidean = blink::raster::create_gdal_raster_from_model
    <int>("output_squared_euclidean.tif", input);
  auto output_manhattan = blink::raster::create_gdal_raster_from_model<int>(
    "output_manhattan.tif", input);
  auto output_chessboard = blink::raster::create_gdal_raster_from_model<int>(
    "output_chessboard.tif", input);

  blink::raster_tools::euclidean_distance_transform(input, output_euclidean, 1);
  blink::raster_tools::chessboard_distance_transform(input, output_chessboard
    , 1);
  blink::raster_tools::manhattan_distance_transform(input, output_manhattan, 1);
  blink::raster_tools::squared_euclidean_distance_transform(input, 
    output_squared_euclidean, 1);
  return 0;
}

int demo_fuzzy_kappa()
{
  auto map1 = blink::raster::open_gdal_raster<int>("map1.rst", GA_ReadOnly);
  auto map2 = blink::raster::open_gdal_raster<int>("map3.rst", GA_ReadOnly);

  auto out = blink::raster::create_gdal_raster_from_model<double>("fk.tif",
    map1);


  //auto mask = open_gdal_raster<int>("region.rst", GA_ReadOnly);
  // create a map of ones, as the mask
  auto nomask = blink::raster::create_temp_gdal_raster_from_model<int>(map1);
  for (auto&& i : nomask)  i = 1;
  const int number_of_categories = 4;
  blink::raster_tools::matrix<double> m(number_of_categories, 
    std::vector<double>(number_of_categories, 0));;
  for (int i = 0; i < number_of_categories; ++i) {
    m[i][i] = 1;
  }
   
  double fuzzykappa;
  bool success = blink::raster_tools::fuzzy_kappa_2009(
    map1,    // input: first map
    map2,    // input: second map
    nomask,   // input: mask map
    number_of_categories, number_of_categories,     // dimension: number of categories in legends
    m,    // parameter: categorical similarity matrix
    blink::raster_tools::exponential_decay{ 2.0 },    // parameter: distance decay function
    out, // result: similarity map
    blink::raster_tools::gdal_raster_maker{}, // RasterMaker::raster<T> r = maker.create<T>(model)
    fuzzykappa);

  std::cout << "Fuzzy Kappa " << fuzzykappa << std::endl;
  return 0;
}


int main()
{
  std::clock_t start = clock();
  
  demo_fuzzy_kappa();
  
  std::clock_t end = clock();
  std::cout << "That took: " << static_cast<double>(end - start) / 
    CLOCKS_PER_SEC << " s" << std::endl;
  return 0;
}