/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

#ifndef SIRF_ABSTRACT_DATA_CONTAINER_TYPE
#define SIRF_ABSTRACT_DATA_CONTAINER_TYPE

#include <map>
#include "sirf/iUtilities/DataHandle.h"

/*
\ingroup Data Container
\brief Abstract data container.

A class for a set of generally heterogeneous items of data.

Has vector features: norm, dot product, linear combination,
which rely on the same features of the items.
*/

namespace sirf {

	typedef std::map<std::string, int> Dimensions;

	class DataContainer {
	public:
		virtual ~DataContainer() {}
		//virtual DataContainer* new_data_container() const = 0;
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const = 0;
		virtual unsigned int items() const = 0;
		virtual float norm() const = 0;
		virtual void dot(const DataContainer& dc, void* ptr) const = 0;
		virtual void multiply
		(const DataContainer& x, const DataContainer& y) = 0;
		virtual void divide
		(const DataContainer& x, const DataContainer& y) = 0;
		virtual void axpby(
			const void* ptr_a, const DataContainer& x,
			const void* ptr_b, const DataContainer& y) = 0;
		bool is_empty() const
		{
			return items() < 1;
		}
		virtual void write(const std::string &filename) const = 0;
		std::unique_ptr<DataContainer> clone() const
		{
			return std::unique_ptr<DataContainer>(this->clone_impl());
		}
	protected:
		virtual DataContainer* clone_impl() const = 0;
	};
}

#endif
