fd = open('_pykdtree.template.cc', 'r')
template = fd.read()
fd.close()

dimensions = range(2,11)
declarations = []

sections = []

dimension_type = dict(name="Double", id="d", type="double")
index_type = dict(id="L", type="unsigned long long")


header = """
//
//  THIS IS A GENERATED FILE DO NOT EDIT
//
//  Edit the gen-pykdtree.py and _pykdtree.template and re-generate
//
//


#include <Python.h>
#include "structmember.h"
#include <kdtree++/kdtree.hpp>
#include <sstream>
#include <iostream>
#include <vector>
#include <limits>

template <size_t DIM, typename COORD_T, typename DATA_T > 
struct record_t {
  static const size_t dim = DIM;
  typedef COORD_T coord_t;
  typedef DATA_T data_t;

  typedef coord_t point_t[dim];
 
  inline coord_t operator[](size_t const N) const { return point[N]; }

  point_t point;
  data_t data;
};


typedef double RANGE_T;


template <class RECORD_T>
inline double tac(RECORD_T r, int k) { return r[k]; }


template <size_t DIM, typename COORD_T, typename DATA_T > 
class PyKDTree {
public:

  typedef record_t<DIM, COORD_T, DATA_T> RECORD_T;
  typedef KDTree::KDTree<DIM, RECORD_T, std::pointer_to_binary_function<RECORD_T,int,double> > TREE_T;
  TREE_T tree;

  PyKDTree() : tree(std::ptr_fun(tac<RECORD_T>)) {  };

  void add(RECORD_T T) { tree.insert(T); };

  /**
     Exact erase.
  */
  bool remove(RECORD_T T) { 
    bool removed = false;

    typename TREE_T::const_iterator it = tree.find_exact(T);
    if (it!=tree.end()) {
      tree.erase_exact(T); 
      removed = true;
    }
    return removed;
  };

  int size(void) { return tree.size(); }

  void optimize(void) { tree.optimise(); }
  
  RECORD_T* find_exact(RECORD_T T) {
    RECORD_T* found = NULL;
    typename TREE_T::const_iterator it = tree.find_exact(T);
    if (it!=tree.end())
      found = new RECORD_T(*it);

    return found;
  }

  size_t count_within_range(typename RECORD_T::point_t T, RANGE_T range) {
    RECORD_T query_record;
    memcpy(query_record.point, T, sizeof(COORD_T)*DIM);

    return tree.count_within_range(query_record, range);
  }

  std::vector<RECORD_T >* find_within_range(typename RECORD_T::point_t T, RANGE_T range) {
    RECORD_T query_record;
    memcpy(query_record.point, T, sizeof(COORD_T)*DIM);

    std::vector<RECORD_T> *v = new std::vector<RECORD_T>;
    tree.find_within_range(query_record, range, std::back_inserter(*v));
    return v;
  }

  RECORD_T* find_nearest (typename RECORD_T::point_t T) {
    RECORD_T* found = NULL;
    RECORD_T query_record;
    memcpy(query_record.point, T, sizeof(COORD_T)*DIM);

    std::pair<typename TREE_T::const_iterator, typename TREE_T::distance_type> best = 
      tree.find_nearest(query_record, std::numeric_limits<typename TREE_T::distance_type>::max());

    if (best.first!=tree.end()) {
      found = new RECORD_T(*best.first);
    }
    return found;
  }

  std::vector<RECORD_T >* get_all() {
    std::vector<RECORD_T>* v = new std::vector<RECORD_T>;

    for (typename TREE_T::const_iterator iter=tree.begin(); iter!=tree.end(); ++iter) {
      v->push_back(*iter);
    }

    return v;
  }

  size_t __len__() { return tree.size(); }
};

"""

upper_defines_3dull = """

#define RECORD_3dull record_t<3, double, unsigned long long >
#define RECORD_3dull_point record_t<3, double, unsigned long long >::point_t
#define RECORD_3dull_coord record_t<3, double, unsigned long long >::coord_t
#define KDTREE_TYPE_3dull KDTree::KDTree<3, RECORD_3dull, std::pointer_to_binary_function<RECORD_3dull,int,double> >

"""

upper_defines = """

#define RECORD_%(NAME)s record_t<%(II)s, %(DIMTYPE)s, %(INDEXTYPE)s >
#define RECORD_%(NAME)s_point record_t<%(II)s, %(DIMTYPE)s, %(INDEXTYPE)s >::point_t
#define RECORD_%(NAME)s_coord record_t<%(II)s, %(DIMTYPE)s, %(INDEXTYPE)s >::coord_t
#define KDTREE_TYPE_%(NAME)s KDTree::KDTree<%(II)s, RECORD_%(NAME)s, std::pointer_to_binary_function<RECORD_%(NAME)s,int,%(DIMTYPE)s> >
#define PYKDTREE_%(NAME)s PyKDTree<%(II)s, %(DIMTYPE)s, %(INDEXTYPE)s >


"""

operators_3dull = """

inline bool operator==(RECORD_3dull const& A, RECORD_3dull const& B) {
    return A.point[0] == B.point[0] && A.point[1] == B.point[1] && A.point[2] == B.point[2] && A.data == B.data;
}

std::ostream& operator<<(std::ostream& out, RECORD_3dull const& T)
{
    return out << '(' << T.point[0] << ',' << T.point[1] << ',' << T.point[2] << '|' << T.data << ')';
}

"""

operators = """

inline bool operator==(%(RECORDNAME)s const& A, %(RECORDNAME)s const& B) {
    return A.data == B.data && %(OPERATOR1)s ;
}

std::ostream& operator<<(std::ostream& out, %(RECORDNAME)s const& T)
{
    return out << '(' %(OPERATOR2)s '|' << T.data << ')';
}


"""

operator1 = "A.point[%i] == B.point[%i]"
operator2 = "<< T.point[%i] <<"


footer_3dull = """

// Global module methods for _pycrates are inserted here
static PyMethodDef module_methods[] = {
  
  //{ (char*) "foo",(PyCFunction)foo, METH_VARARGS, (char*) "foo doc string"},
  
  {NULL, NULL, 0, NULL}  // Sentinel
};


#ifndef PyMODINIT_FUNC	// declarations for DLL import/export
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_pykdtree(void)
{ 
  PyObject *m;

  if (PyType_Ready(&PyTree_3dullType) < 0)
    return;

  //import_array();

  m = Py_InitModule3((char*)"_pykdtree", module_methods,
		     (char*)"Example module that creates an extension type.");
  
  if (m == NULL)
    return;
  
  Py_INCREF(&PyTree_3dullType);

  // C++ classes mapped to Python types go here and are then inserted into
  // the module namespace as Python classes
  PyModule_AddObject(m, (char*)"PyTree_3dull", (PyObject *)&PyTree_3dullType);
}


"""

footer = """

// Global module methods for _pycrates are inserted here
static PyMethodDef module_methods[] = {
  
  //{ (char*) "foo",(PyCFunction)foo, METH_VARARGS, (char*) "foo doc string"},
  
  {NULL, NULL, 0, NULL}  // Sentinel
};


#ifndef PyMODINIT_FUNC	// declarations for DLL import/export
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_pykdtree(void)
{ 
  PyObject *m;

  %(READYNAME)s

  //import_array();

  m = Py_InitModule3((char*)"_pykdtree", module_methods,
		     (char*)"Example module that creates an extension type.");
  
  if (m == NULL)
    return;
  
  %(INCREFNAME)s

  // C++ classes mapped to Python types go here and are then inserted into
  // the module namespace as Python classes
  %(ADDNAME)s
}


"""



pages = []
readynames = []
increfnames = []
addnames = []

for ii in dimensions:

    name = "%i%s%s" % (ii, dimension_type["id"], index_type["id"])
    dimids = "".join([dimension_type["id"] for jj in range(ii)])





    udefines = str(upper_defines)
    udefines = udefines.replace("%(NAME)s", name)
    udefines = udefines.replace("%(II)s", str(ii))
    udefines = udefines.replace("%(DIMTYPE)s", dimension_type["type"])
    udefines = udefines.replace("%(INDEXTYPE)s", index_type["type"])

    structname = "PyTree_%s" % name
    treename = "PYKDTREE_%s" % name
    recordname = "RECORD_%s" % name
    pointname = "RECORD_%s_point" % name
    coordname = "RECORD_%s_coord" % name

    operator1s = []
    operator2s = []

    iter_points = []
    temp_points = []
    points = []
    record_points = []



    for jj in range(ii):
        
        iter_points.append("(*iter).point[%i]" % jj)
        temp_points.append("temp->point[%i]" % jj)
        points.append("&point_t[%i]," % jj)
        record_points.append("&record.point[%i]," % jj)

        operator1s.append(operator1 % (jj, jj))
        operator2s.append(operator2 % (jj))

    operator1s = " && ".join(operator1s)
    operator2s = " ',' ".join(operator2s)

    operator = str(operators)
    operator = operator.replace("%(RECORDNAME)s", recordname)
    operator = operator.replace("%(OPERATOR1)s", operator1s)
    operator = operator.replace("%(OPERATOR2)s", operator2s)

    iter_points = ",".join(iter_points)
    temp_points = ",".join(temp_points)
    points = "\n".join(points)
    record_points = "\n".join(record_points)

    classname = "KDTree_%i%s" % (ii, dimension_type["name"])

    page = str(template)
    page = page.replace("%(STRUCTNAME)s", structname)
    page = page.replace("%(CLASSNAME)s", classname)
    page = page.replace("%(TREENAME)s", treename)
    page = page.replace("%(RECORDNAME)s", recordname)
    page = page.replace("%(POINTNAME)s", pointname)
    page = page.replace("%(COORDNAME)s", coordname)
    page = page.replace("%(DIMIDS)s", dimids)
    page = page.replace("%(INDEXID)s", index_type["id"])
    page = page.replace("%(II)s", str(ii))
    page = page.replace("%(DIMTYPE)s", dimension_type["type"])
    page = page.replace("%(INDEXTYPE)s", index_type["type"])
    page = page.replace("%(ITERPOINTS)s", iter_points)
    page = page.replace("%(TEMPPOINTS)s", temp_points)
    page = page.replace("%(RECORDPOINTS)s", record_points)
    page = page.replace("%(POINTS)s", points)
    page = page.replace("%(POINTS1)s", points.strip(","))

    pages.append(udefines + operator + page)

    readynames.append("""if (PyType_Ready(&%sType) < 0)
    return;
""" % structname)
    increfnames.append("Py_INCREF(&%sType);" % structname)
    addnames.append('PyModule_AddObject(m, (char*)"%s", (PyObject *)&%sType);' % (classname, structname))



footer = footer.replace("%(READYNAME)s", "\nelse ".join(readynames))
footer = footer.replace("%(INCREFNAME)s", "\n".join(increfnames))
footer = footer.replace("%(ADDNAME)s", "\n".join(addnames))

page = header + "".join(pages) + footer
fd = file("_pykdtree.cc", "w")
fd.write(page)
fd.close()
