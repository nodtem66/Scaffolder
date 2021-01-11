#define SOL_ALL_SAFETIES_ON 1
#include "implicit_function.h"
#include <sol/sol.hpp>


int main(int argc, char* argv[]) {
	std::cout << "=== opening a state ===" << std::endl;
	std::cout << argc << std::endl;

	if (argc >= 2) {
		sol::state lua;
		// open some common libraries
		lua.open_libraries(sol::lib::base, sol::lib::math, sol::lib::string);
		sol::load_result lua_file = lua.load_file(argv[1]);
		if (!lua_file.valid()) {
			sol::error err = lua_file;
			std::cerr << "[Lua Error] " << err.what() << std::endl;
		}
		else {
			lua.set("w", 1);
			lua.set("t", 2);
			lua.set("grid_size", 3);
			lua.set_function("v", [](double x, double y, double z) -> double {return 2 * x; });
			lua_file();
			sol::function fn = lua["surface"];
			Function* surface;
			surface = new LuaFunction(fn);
			Implicit_function f(surface, 1);
			std::cout << (double)(f)(10, 10, 10) << std::endl;
		}
	}
	

	return 0;
}