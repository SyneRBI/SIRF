//
//
/*!
  \file 
 
  \brief defines the stir::error() function

  \author Kris Thielemans
  \author PARAPET project


*/
/*
    Copyright (C) 2000 PARAPET partners
    Copyright (C) 2000- 2010, Hammersmith Imanet Ltd
    This file is part of STIR.

    This file is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    See STIR/LICENSE.txt for details
*/
#include "stir/error.h"

/* EO 22/07/15 */
#include "TextWriter.h"
#include "StirException.h"
/* EO 25/08/15 */
#define cerr Stir::cerr
#define endl Stir::endl

#include <cstdarg>
#include <string>
#include <iostream>
/* Warning: vsnprintf is only ISO C99. So your compiler might not have it.
   Visual Studio can be accomodated with the following work-around
*/
#ifdef BOOST_MSVC
#define vsnprintf _vsnprintf
#endif

START_NAMESPACE_STIR

void error(const char *const s, ...)
{  
  va_list ap;
  va_start(ap, s);
  const unsigned size=10000;
  char tmp[size];
  const int returned_size= vsnprintf(tmp,size, s, ap);
  if (returned_size<0)
		cerr << "\nERROR: but error formatting error message" << endl;
		//std::cerr << "\nERROR: but error formatting error message" << std::endl;
	else
    {
			///* EO 22/07/15 */
			//writeText("\nERROR: ");
			//writeText(tmp);
			//writeText("\n");
			cerr << "\nERROR: " << tmp << endl;
			//std::cerr << "\nERROR: " << tmp <<std::endl;
			if (static_cast<unsigned>(returned_size) >= size)
				cerr << "\nWARNING: previous error message truncated as it exceeds "
				<< (int)size << "bytes" << endl;
			//std::cerr << "\nWARNING: previous error message truncated as it exceeds "
			//	<< size << "bytes" << std::endl;
	}
  std::string msg = tmp;
  //throw msg;
	throw StirException(tmp, __FILE__, __LINE__);
}
END_NAMESPACE_STIR
