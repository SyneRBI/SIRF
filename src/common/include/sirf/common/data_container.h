/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

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

#ifndef SIRF_DATA_CONTAINER_TYPE
#define SIRF_DATA_CONTAINER_TYPE

/*
\ingroup Data Container
\brief Abstract data container.

A class for a set of generally heterogeneous items of data.

Has vector features: norm, dot product, linear combination,
which rely on the same features of the items.
*/
template <typename T>
class aDataContainer {
public:
	virtual ~aDataContainer() {}
	virtual aDataContainer<T>* new_data_container() = 0;
	virtual unsigned int items() = 0;
	virtual float norm() = 0;
	virtual T dot(const aDataContainer<T>& dc) = 0;
	virtual void multiply
		(const aDataContainer<T>& x, const aDataContainer<T>& y) = 0;
	virtual void divide
		(const aDataContainer<T>& x, const aDataContainer<T>& y) = 0;
	virtual void axpby(
		T a, const aDataContainer<T>& x,
		T b, const aDataContainer<T>& y) = 0;
};

#endif