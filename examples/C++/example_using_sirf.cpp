#include <iostream>
#include "sirf/STIR/stir_data_containers.h"

int main(int argc, char**argv)
{
  if (argc != 2)
    {
      std::cerr << "Need an image filename as argument\n";
      return EXIT_FAILURE;
    }
  std::string filename(argv[1]);
  sirf::STIRImageData im(filename);
  std::cout <<"geometrical info:\n" << im.get_geom_info_sptr()->get_info();

  return EXIT_SUCCESS;
}
