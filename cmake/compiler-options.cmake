include(apple-llvm-toolchain)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") # Matches "Clang" and "AppleClang"
	if (APPLE)
		message(STATUS "Using clang ${CMAKE_CXX_COMPILER_VERSION} on apple")
		set(CMAKE_MACOSX_RPATH 1)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -stdlib=libc++ -pthread")
	else()
		if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
			message(FATAL_ERROR " clang version must be at least 5!")
		endif()
		message(STATUS "Using clang ${CMAKE_CXX_COMPILER_VERSION}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread")
	endif()
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fomit-frame-pointer -msse -msse2 -funroll-loops")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O1")
	if (NOT "${CLANG_SANITIZER}" STREQUAL "none")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=${CLANG_SANITIZER}")
	endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7)
		message(FATAL_ERROR " g++ version must be at least 7!")
	endif()
	message(STATUS "Using g++ ${CMAKE_CXX_COMPILER_VERSION}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunknown-pragmas")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread -fmax-errors=5 -fdiagnostics-color=auto")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O1")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	message("-- Using Visual Studio C++ ${CMAKE_CXX_COMPILER_VERSION}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP1 /D_SCL_SECURE_NO_WARNINGS /std:c++17")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /MT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd /bigobj")
else()
	message(WARNING "Possibly unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

if(DEVELOPER)
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wconversion")
	endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG -ftemplate-backtrace-limit=0")
	#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -print-search-dirs")

	if (ALLWARNINGS)
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++98-compat -Wno-c++98-compat-pedantic")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-old-style-cast")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reserved-id-macro")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-newline-eof")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-documentation")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-documentation-unknown-command")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-weak-vtables")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-global-constructors")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-exit-time-destructors")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-switch-enum")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-covered-switch-default")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-padded")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-float-equal")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedef")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-variable-declarations")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory")
	endif()
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
	endif()
	endif()

	set(CMAKE_BUILD_TYPE "DEBUG")
endif()

message(STATUS "CXX Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "CXX Debug Flags: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CXX Release Flags: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
