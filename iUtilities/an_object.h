#ifndef ROOT_OBJECT_TYPE
#define ROOT_OBJECT_TYPE

#include <string>

class anObject {
public:
	std::string class_name() const
	{
		return class_;
	}
protected:
	std::string class_;
};

#endif
