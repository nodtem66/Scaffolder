#ifndef UTILS_INCLUDED
#define UTILS_INCLUDED
#include <iostream>
#include <string>
#include <locale>

#ifdef _WIN32
#include <Windows.h>
#include <direct.h>
#else
#include <sys/time.h>
#include <ctime>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

namespace util {

    class NullBuffer : public std::streambuf
    {
    public:
        int overflow(int c) { return c; }
    };

    class NullStream : public std::ostream {
    public:
        NullStream() : std::ostream(&m_sb) {}
        static NullStream& getInstance() {
            static NullStream s;
            return s;
        }
    private:
        NullBuffer m_sb;
    };

    typedef long long int64; typedef unsigned long long uint64;

    int make_dir(std::string& str);

    util::uint64 GetTimeMs64();

    std::string PathGetBasename(std::string const& path);

    std::string PathGetExtension(std::string const& path);

    void to_lower(std::string& s);
}

#if defined(__clang__) && (__cplusplus >= 201500L)
#include <algorithm>
#include <random>
namespace std {
    template <class RandomAccessIterator>
    inline void random_shuffle(RandomAccessIterator first, RandomAccessIterator last);
    
    template <class RandomAccessIterator, class RandomNumberGenerator>
    inline void random_shuffle(RandomAccessIterator first, RandomAccessIterator last,
        RandomNumberGenerator& gen);

    template <class RandomAccessIterator, class RandomNumberGenerator>
    inline void random_shuffle(RandomAccessIterator first, RandomAccessIterator last,
        RandomNumberGenerator&& gen);
}
#endif
#endif