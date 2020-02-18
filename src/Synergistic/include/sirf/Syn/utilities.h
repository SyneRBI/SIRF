/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 Rutherford Appleton Laboratory STFC
Copyright 2020 University College London

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

/*!
\file
\ingroup SIRF C++ Utilities
\brief C++ utilities.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#pragma once

namespace sirf {
	class ImageDataWrap {
	public:
		ImageDataWrap(const std::string &filename, const std::string &engine, bool verbose);
		ImageData& data()
		{
			return *img_sptr_;
		}
		const ImageData& data() const
		{
			return *img_sptr_;
		}
		std::shared_ptr<ImageData> data_sptr()
		{
			return img_sptr_;
		}
		const std::shared_ptr<ImageData> data_sptr() const
		{
			return img_sptr_;
		}
	private:
		std::shared_ptr<ImageData> img_sptr_;
	};
}