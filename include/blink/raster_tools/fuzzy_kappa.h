//
//=======================================================================
// Copyright 2016
// Author: Alex Hagen-Zanker
// University of Surrey
//
// Distributed under the MIT Licence (http://opensource.org/licenses/MIT)
//=======================================================================
//
// This implements the Fuzzy Kappa method

#ifndef BLINK_RASTER_TOOLS_FUZZY_KAPPA_H_AHZ
#define BLINK_RASTER_TOOLS_FUZZY_KAPPA_H_AHZ

#include <blink/raster_tools/distance_transform.h>
#include <blink/raster/raster_traits.h>
#include <blink/raster/utility.h>

#include <blink/iterator/range_range.h>
#include <blink/iterator/zip_range.h>

#include <algorithm>
#include <vector> // vectors are used to store arrays
#include <map>    // maps are used to create distributions
namespace blink {
  namespace raster_tools {
    typedef std::map<double, int, std::greater<double> > distribution;

    template<class T>
    using matrix = std::vector< std::vector<T> >;

    ////////////////////////////////////////////////////////////////////////////////
    // Fuzzy Kappa can be based on any distance decay function. The function may 
    // be set by one or more parameters. Therefore it is passed as a functor
    //
    class exponential_decay 
    {
    public:
      exponential_decay(double halving) : m_halving(halving)
      {
      }

      double operator()(double d) const
      {
        return pow(0.5, d / m_halving);
      }

      double m_halving;
    };

    class one_neighbour
    {
    public:
      one_neighbour(double value) : m_value(value)
      {
      }

      double operator()(double d) const
      {
        if (d < 0.9) return 1.0;
        if (d < 1.1) return m_value;
        return 0;
      }
      double m_value;

    };


    class gdal_raster_maker
    {
    public:
      template<class T>
      using raster_type = blink::raster::gdal_raster<T>;

      template<class T, class U>
      raster_type<T> create(const blink::raster::gdal_raster<U>& model)
      {
        return blink::raster::create_temp_gdal_raster_from_model<T>(model);
      }
    };

    ////////////////////////////////////////////////////////////////////////////////
    // This function takes two distribution and returns the expected minimum value 
    // when a number is sampled from both functions
    //
    double expected_minumum_of_two_distributions(
      const distribution& distriA, // pairs of values and counts
      const distribution& distriB,
      double totalA,               //total count
      double totalB)
    {
      //std::map::const_reverse_iterator goes through the elements from low to high
      distribution::const_iterator iterA = distriA.begin();
      distribution::const_iterator iterB = distriB.begin();

      double pCum = 0;         // cumulative probability p(min(A,B) <= x)
      double sumA = 0;         // intermediate for pCumA
      double sumB = 0;         // intermediate for pCumB
      double pCumA = 0;        // cumulative probabilty p(A<=x)
      double pCumB = 0;        // cumulative probability p(B <=x)
      double x = 1;            // x
      double expected = 0;     //Integral of x * p(x) over x

      while (iterA != distriA.end() || iterB != distriB.end()) {
        // In order to pass all values in both distributions advance to the 
        // largest next value, except if already at the end of the series 
        if (iterB == distriB.end()
          || (iterA != distriA.end() && iterA->first > iterB->first)) {
          x = iterA->first;
          sumA += iterA->second;
          pCumA = sumA / totalA;
          iterA++;
        }
        else {
          x = iterB->first;
          sumB += iterB->second;
          pCumB = sumB / totalB;
          iterB++;
        }

        const double pCumPrevious = pCum;
        pCum = pCumA * pCumB; // p(A>=x AND B>=x )
        expected += (pCum - pCumPrevious) * x;
      }
      return expected;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // This is the entry function to calculate Fuzzy Kappa (improved), 
    // Returns false if there are no cells to compare (Fuzzy Kappa = 0)
    // Returns false if all cells in both maps are identical (Fuzzy Kappa = 1)
    //
    template<class RasterA, class RasterB, class RasterMask, class RasterOut,
    class DistanceDecay, class RasterMaker>
      bool fuzzy_kappa_2009(
      RasterA& mapA,              // input: first map
      RasterB& mapB,              // input: second map
      RasterMask& mask,           // input: mask map
      int nCatsA, int nCatsB,     // dimension: number of categories in legends
      const matrix<double>& m,    // parameter: categorical similarity matrix
      DistanceDecay f,            // parameter: distance decay function
      RasterOut& comparison,      // result: similarity map
      RasterMaker maker,          // RasterMaker::raster<T> r = maker.create<T>(model)
      double& fuzzykappa)         // result: improved fuzzy kappa
    {
        // for use with std::get
        using a_value = blink::raster::raster_traits::value_type<RasterA>;
        using b_value = blink::raster::raster_traits::value_type<RasterB>;
        using mask_value = blink::raster::raster_traits::value_type<RasterMask>;
        using out_value = blink::raster::raster_traits::value_type<RasterOut>;
        using temp_raster = RasterMaker::raster_type<double>;

        // Nearest neighbour distances for all categories in both maps
        std::vector<temp_raster> distancesA;
        std::vector<temp_raster> distancesB;
        for (int catA = 0; catA < nCatsA; ++catA) {
          temp_raster temp = maker.create<double>(mapA);
          blink::raster_tools::euclidean_distance_transform(mapA, temp, catA);
          //Apply distance decay function on the distances
          for (auto&& i : temp) i = f(i);
          distancesA.emplace_back(std::move(temp));
        }

        for (int catB = 0; catB < nCatsB; catB++) {
          temp_raster temp = maker.create<double>(mapB);
          blink::raster_tools::euclidean_distance_transform(mapB, temp, catB);
          //Apply distance decay function on the distances
          for (auto&& i : temp) i = f(i);
          distancesB.emplace_back(std::move(temp));
        }

        // Apply the categorical similarity matrix. Note that this yields the 
        // similarity of map A to the categories of map B and vice versa
        std::vector<temp_raster> similarityA;
        std::vector<temp_raster> similarityB;
        for (int catA = 0; catA < nCatsA; ++catA) {
          auto temp = maker.create<double>(mapB);
          for (auto&& i : temp) i = 0;
          similarityB.emplace_back(std::move(temp));
        }
        for (int catB = 0; catB < nCatsB; ++catB) {
          auto temp = maker.create<double>(mapB);
          for (auto&& i : temp) i = 0;
          similarityA.emplace_back(std::move(temp));
        }

        auto multi_similarity_a = blink::iterator::make_range_range(similarityA);
        auto multi_similarity_b = blink::iterator::make_range_range(similarityB);
        auto multi_distance_a = blink::iterator::make_range_range(distancesA);
        auto multi_distance_b = blink::iterator::make_range_range(distancesB);
        static const int a_index = 0;
        static const int b_index = 1;
        static const int mask_index = 2;
        static const int out_index = 3;
        static const int sim_a_index = 4;
        static const int sim_b_index = 5;
        static const int dist_a_index = 6;
        static const int dist_b_index = 7;

        auto zip = blink::iterator::make_zip_range(
          std::ref(mapA), //0
          std::ref(mapB), //1
          std::ref(mask), //2
          std::ref(comparison), //3
          std::ref(multi_similarity_a), //4
          std::ref(multi_similarity_b), //5
          std::ref(multi_distance_a),//6
          std::ref(multi_distance_b)//7
          );

        double mean = 0;
        int count = 0;
        std::vector< std::vector<distribution> > distributionA(nCatsA,
          std::vector<distribution>(nCatsB));
        std::vector< std::vector<distribution> > distributionB(nCatsB,
          std::vector<distribution>(nCatsA));

        // count number of cell and cells per category in both maps
        std::vector<int> catCountsA(nCatsA, 0);
        std::vector<int> catCountsB(nCatsB, 0);

        for (auto&& i : zip) {
          for (int catA = 0; catA < nCatsA; catA++){
            for (int catB = 0; catB < nCatsB; catB++){
              const double mAB = m[catA][catB];

              const double da = std::get<dist_a_index>(i)[catA];
              const double db = std::get<dist_b_index>(i)[catB];
              const double simAB = mAB * da;
              const double simBA = mAB * db;

              const double sa = std::get<sim_a_index>(i)[catB];
              const double sb = std::get<sim_b_index>(i)[catA];

              if (sa < simAB) std::get<sim_a_index>(i)[catB] = simAB;
              if (sb < simBA) std::get<sim_b_index>(i)[catA] = simBA;
            }
          }
          if (std::get<2>(i)) {
            const a_value catA = std::get<a_index>(i);
            const b_value catB = std::get<b_index>(i);
            ++catCountsA[catA];
            ++catCountsB[catA];

            const double simAB = std::get<sim_a_index>(i)[catB];
            const double simBA = std::get<sim_b_index>(i)[catA];
            const double sim = std::min(simAB, simBA);
            mean += sim;
            ++count;
            std::get<out_index>(i) = sim;

            for (int catB = 0; catB < nCatsB; catB++) {
              const double sim = std::get<sim_a_index>(i)[catB];
              distributionA[catA][catB][sim]++; // See documentation of std::map
            }

            for (int catA = 0; catA < nCatsA; catA++) {
              const double sim = std::get<sim_b_index>(i)[catA];
              distributionB[catB][catA][sim]++; // See documentation of std::map
            }
          }
          else {
            std::get<out_index>(i) = -1; //nodata value
          }
        }

        if (count == 0) {
          fuzzykappa = 0;
          return false;
        }
        mean /= count;

        // Calculate expected similarity
        double expected = 0;
        const double squaredTotal = (double)(count)*(double)(count);
        for (int catA = 0; catA < nCatsA; catA++){
          for (int catB = 0; catB < nCatsB; catB++){
            // The if statement avoids division by zero
            if (catCountsA[catA] > 0 && catCountsB[catB] > 0) {
              const double pCats = (double)(catCountsA[catA]) * (double)(catCountsB[catB])
                / squaredTotal;
              const double eCats = expected_minumum_of_two_distributions(
                distributionA[catA][catB], distributionB[catB][catA], catCountsA[catA],
                catCountsB[catB]);

              expected += pCats * eCats;
            }
          }
        }

        // If all cells are identical to each other
        if (expected == 1) {
          fuzzykappa = 1;
          return false;
        }

        // Calculate Fuzzy Kappa 
        fuzzykappa = (mean - expected) / (1.0 - expected);
        //std::cout << "Mean " << mean << std::endl;
        //std::cout << "Expected " << expected << std::endl;
        return true;
    }
  }
}
#endif
