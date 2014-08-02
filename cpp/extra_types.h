/// \file extra_types.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// Provides some extra types that may be generally useful

#if !defined(_XDRESS_EXTRA_TYPES_)
#define _XDRESS_EXTRA_TYPES_

#if defined(__cplusplus)
namespace extra_types
{
  /// complex type struct, matching PyTables definition
//  typedef struct {
//    double re;  ///< real part
//    double im;  ///< imaginary part
//  } complex_t;

  /// Chivalrously handles C++ memory issues that Cython does
  /// not yet have a syntax for.  This is a template class,
  /// rather than three template functions, because Cython does
  /// not yet support template function wrapping.
  template <class T>
  class MemoryKnight
  {
    public:
      MemoryKnight(){};   ///< Default constructor
      ~MemoryKnight(){};  ///< Default Destructor

      /// Creates a new instance of type T on the heap using
      /// its default constructor.
      /// \return T *
      T * defnew(){return new T();};

      /// Creates a new instance of type T, using T's default 
      /// constructor, at a given location.
      /// \param void * ptr, location to create T instance
      /// \return value of ptr recast as T *
      T * renew(void * ptr){return new (ptr) T();};

      /// Deallocates a location in memory using delete. 
      /// \param T * ptr, location to remove
      void deall(T * ptr){delete ptr;};
  };

// End namespace extra_types
};

#elif defined(__STDC__)

// de nada

#endif


/// complex type struct, matching PyTables definition
typedef struct {
  double re;  ///< real part
  double im;  ///< imaginary part
} xd_complex_t;

#endif

