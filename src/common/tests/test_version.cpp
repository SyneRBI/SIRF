#include "version.h"
#include <iostream>


int main(int argc, char ** argv){

    std::cout << SIRF_VERSION_MAJOR << "." << SIRF_VERSION_MINOR << std::endl;
    
    return 0;

}