//
//=======================================================================
// Copyright 2016
// Author: Alex Hagen-Zanker
// University of Surrey
//
// Distributed under the MIT Licence (http://opensource.org/licenses/MIT)
//=======================================================================
//
// This implements the distance transform method by Meijster. 
// The method is very amenable to parallelization. 
// But not done here..

#ifndef BLINK_RASTER_TOOLS_DISTANCE_TRANSFORM_H_AHZ
#define BLINK_RASTER_TOOLS_DISTANCE_TRANSFORM_H_AHZ


#include <blink/raster/raster_traits.h>
#include <algorithm>
#include <vector>




 namespace blink {
   namespace raster_tools {
    namespace detail
    {
      struct euclidean{};
    }
    struct euclidean_squared : public detail::euclidean{};
    struct euclidean_non_squared : public detail::euclidean{};
    struct manhattan{};
    struct chessboard{};

    namespace detail
    {
       int round(const double& f)
      {
        return static_cast<int>(f + 0.5);
      }

      template<class T>
      double optionally_square_root(const T& value, 
        const euclidean_non_squared&)
      {
        return sqrt(static_cast<double>(value));
      }

      template<class T, class OtherMethod>
      T optionally_square_root(const T& value, const OtherMethod&)
      {
        return value;
      }
      int f(int x, int i, std::vector<int>& g, const detail::euclidean&)
      {
        const int dx = x - i;
        const int dy = g[i];
        return dx * dx + dy * dy;
      }

      int f(int x, int i, std::vector<int>& g, const manhattan&)
      {
        return abs(x - i) + g[i];
      }

      int f(int x, int i, std::vector<int>& g, const chessboard&)
      {
        return std::max(abs(x - i), g[i]);
      }

      int sep(int i, int u, std::vector<int>& g, int, const detail::euclidean&)
      {
        return ( (u-i) * (u+i) + (g[u] - g[i]) * (g[u] + g[i] ) ) / (2 * (u - i));

        //return (u*u - i*i + g[u] * g[u] - g[i] * g[i]) / (2 * (u - i));
      }

      int sep(int i, int u, std::vector<int>& g, int inf, const manhattan&)
      {
        if (g[u] >= g[i] + u - i) return inf;
        if (g[i] > g[u] + u - i) return -inf;
        return (g[u] - g[i] + u + i) / 2;
      }

      int sep(int i, int u, std::vector<int>& g, int, const chessboard&)
      {
        if (g[i] <= g[u]) return std::max(i + g[u], (i + u) / 2);
        return std::min(u - g[i], (i + u) / 2);
      }

      struct st_pair
      {
        st_pair(int s, int t) : s(s), t(t)
        {}
        int s;
        int t;
      };

      template<class T, class ResultIter, class MethodTag>
      void process_line(std::vector<T>& g, ResultIter& iter, int inf, 
        const MethodTag&)
      {
        const int m = g.size();
        std::vector<st_pair> st(1, st_pair(0, 0));
        for (int u = 1; u < m; ++u) {
          while (!st.empty() &&
            f(st.back().t, st.back().s, g, MethodTag{})
        > f(st.back().t, u, g, MethodTag{})){
            st.pop_back();
          }
          if (st.empty()){
            st.emplace_back(u, 0);
          }
          else {
            const int w = 1 + sep(st.back().s, u, g, inf, MethodTag{});
            if (w < m){
              st.emplace_back(u, w);
            }
          }
        }
        for (int u = m - 1; u >= 0; --u, ++iter) { 
          //++iter because g was in reverse
          *iter = optionally_square_root(f(u, st.back().s, g, MethodTag{}), 
            MethodTag{});
          if (u == st.back().t){
            st.pop_back();
          }
        }
      }
    }
    
    template<class InRaster, class OutRaster>
    void euclidean_distance_transform(const InRaster& in, OutRaster& out, 
      const blink::raster::raster_traits::value_type<InRaster>& target)
    {
      distance_transform(in, out, target, euclidean_non_squared{});
    }
    
    template<class InRaster, class OutRaster>
    void squared_euclidean_distance_transform(const InRaster& in, OutRaster& out
      , const blink::raster::raster_traits::value_type<InRaster>& target)
    {
      distance_transform(in, out, target, euclidean_squared{});
    }

    template<class InRaster, class OutRaster>
    void manhattan_distance_transform(const InRaster& in, OutRaster& out, 
      const blink::raster::raster_traits::value_type<InRaster>& target)
    {
      distance_transform(in, out, target, manhattan{});
    }

    template<class InRaster, class OutRaster>
    void chessboard_distance_transform(const InRaster& in, OutRaster& out,
      const blink::raster::raster_traits::value_type<InRaster>& target)
    {
      distance_transform(in, out, target, chessboard{});
    }

    template<class InRaster, class OutRaster, class Method>
    void distance_transform(const InRaster& in, OutRaster& out,
      const blink::raster::raster_traits::value_type<InRaster>& target,
      const Method&)
    {
      using in_type = blink::raster::raster_traits::value_type<InRaster>;
      using out_type = blink::raster::raster_traits::value_type<OutRaster>;
      const int rows = blink::raster::raster_operations::size1(in);
      const int cols = blink::raster::raster_operations::size2(in);

      const out_type inf = rows + cols;
      
      auto a = in.begin();
      auto v = out.begin();

      // Going by rows instead of by columns. This should be changed in a 
      // parallel implementation
      
      for (int c = 0; c < cols; ++c, ++a, ++v) { // first row
        (*v) = (static_cast<in_type>(*a) == target) ? 0 : inf;
      }

      auto u = out.begin();
      for (int r = 1; r < rows; ++r) {
        for (int c = 0; c < cols; ++c, ++a, ++u, ++v) { // subsequent rows
          out_type up = *u;
          (*v) = (static_cast<in_type>(*a) == target) ? 0 : up == inf ? inf : up + 1;
        }
      }

      u = out.begin() + ((rows - 1) * cols - 1);
      v = out.begin() + (rows * cols - 1);

      for (int r = rows - 2; r >= 0; --r) {
        std::vector<int> g;
        for (int c = 0; c < cols; ++c, --u, --v) { // going back
          out_type up = (*u);
          out_type vp = (*v);
          if (up > vp) {
            (*u) = vp + 1;
          }
          g.push_back(round(vp)); // will have values conveniently in reverse
        }
        detail::process_line(g, out.begin() + (r + 1) * cols, inf, Method{});
      }
      std::vector<int> g; //last line
      for (int c = 0; c < cols; ++c, --v) { // going back
        out_type vp = (*v);
        g.push_back(round(vp)); // will have values conveniently in reverse
      }
      detail::process_line(g, out.begin(), inf, Method{});
    }
  }
}
#endif