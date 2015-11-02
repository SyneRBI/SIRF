#ifndef TRY_STATIC
#define TRY_STATIC

#include <iostream>

class A {
public:
	static void init_i() {
		static bool initialized = false;
		if (!initialized) {
			i_ = -100;
			initialized = true;
		}
	}
	void set_i(int i) {
		init_i();
		i_ = i;
	}
	int get_i() const {
		return i_;
	}
	void print_i() {
		init_i();
		std::cout << i_ << std::endl;
	}
private:
	static int i_;
};

#endif
