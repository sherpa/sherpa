//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_
#ifndef __sherpa_extension_hh__
#define __sherpa_extension_hh__

#include <sherpa/array.hh>


typedef unsigned int SherpaUInt;
typedef sherpa::Array< unsigned int, NPY_UINT > SherpaUIntArray;

typedef double SherpaFloat;
typedef sherpa::Array< double, NPY_DOUBLE > DoubleArray;
typedef DoubleArray SherpaFloatArray;

typedef sherpa::Array< char, NPY_CHAR > CharArray;
typedef sherpa::Array< int, NPY_INT > IntArray;

typedef int (*converter)( PyObject*, void* );

#define CONVERTME(arg) ((converter) sherpa::convert_to_contig_array<arg>)

#define SHERPAMOD(name, fctlist) \
PyMODINIT_FUNC \
init##name(void) \
{ \
  import_array(); \
  Py_InitModule( (char*)#name, fctlist ); \
}


#define FCTSPEC(name, func) \
 { (char*)#name, (PyCFunction)func, METH_VARARGS, NULL }


#endif /* __sherpa_extension_hh__ */
