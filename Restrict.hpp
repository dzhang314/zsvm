#ifndef ZSVM_RESTRICT_HPP_INCLUDED
#define ZSVM_RESTRICT_HPP_INCLUDED

#if defined(_MSC_VER)
#define RESTRICT __restrict
#elif defined(__GNUC__)
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

#endif // ZSVM_RESTRICT_HPP_INCLUDED
