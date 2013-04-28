#include <iostream>
#include <boost/ptr_container/ptr_vector.hpp>

#include "particles.hpp"

int main()
{
  boost::ptr_vector<particles_proto<float>> ptrs;

  std::cerr << "allocating OpenMP..." << std::endl;
  ptrs.push_back(new particles<float, openmp>());

  std::cerr << "allocating CUDA..." << std::endl;
  ptrs.push_back(new particles<float, cuda>());

  std::cerr << "looping (using C++11)..." << std::endl;
  for (auto &ptr : ptrs)
  {
    ptr.func();
  }

  std::cerr << "done." << std::endl;
}