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

#include <fstream>
#include <string.h>

#include "stir/TextWriter.h"

#include "cstir_shared_ptr.h"
#include "data_handle.h"

using namespace SPTR_NAMESPACE;

USING_NAMESPACE_STIR

extern "C" {
	void* newTextPrinter(const char* stream) 
	{
		TextPrinter* p = new TextPrinter(stream);
		return (void*)p;
	}
	void* newTextWriter(const char* stream) 
	{
		TextWriter* w = new TextWriter;
		if (strlen(stream)) {
			w->out = new std::ofstream;
			((std::ofstream*)w->out)->open(stream, std::ios::out);
		}
		return (void*)w;
	}
	void openChannel(int channel, void* ptr_w) 
	{
		TextWriterHandle h;
		switch (channel) {
		case INFORMATION_CHANNEL:
			h.set_information_channel((aTextWriter*)ptr_w);
			break;
		case WARNING_CHANNEL:
			h.set_warning_channel((aTextWriter*)ptr_w);
			break;
		case ERROR_CHANNEL:
			h.set_error_channel((aTextWriter*)ptr_w);
			break;
		default:
			h.set_information_channel((aTextWriter*)ptr_w);
			h.set_warning_channel((aTextWriter*)ptr_w);
			h.set_error_channel((aTextWriter*)ptr_w);
		}
	}
	void closeChannel(int channel, void* ptr_w) 
	{
		TextWriterHandle h;
		switch (channel) {
		case INFORMATION_CHANNEL:
			if (h.information_channel_ptr() == ptr_w)
				h.set_information_channel(0);
			break;
		case WARNING_CHANNEL:
			if (h.warning_channel_ptr() == ptr_w)
				h.set_warning_channel(0);
			break;
		case ERROR_CHANNEL:
			if (h.error_channel_ptr() == ptr_w)
				h.set_error_channel(0);
			break;
		default:
			if (h.information_channel_ptr() == ptr_w)
				h.set_information_channel(0);
			if (h.warning_channel_ptr() == ptr_w)
				h.set_warning_channel(0);
			if (h.error_channel_ptr() == ptr_w)
				h.set_error_channel(0);
		}
	}
	void* deleteTextPrinter(void* ptr) 
	{
		try {
			delete (TextPrinter*)ptr;
			return new DataHandle;
		}
		CATCH;
	}
	void* deleteTextWriter(void* ptr_w) 
	{
		try {
			if (!ptr_w)
				return new DataHandle;
			TextWriter* w = (TextWriter*)ptr_w;
			if (!w->out)
				return new DataHandle;
			((std::ofstream*)w->out)->close();
			((std::ofstream*)w->out)->clear();
			delete w->out;
			delete w;
			return new DataHandle;
		}
		CATCH;
	}
}
