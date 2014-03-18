
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



#define RECORD_2dL record_t<2, double, unsigned long long >
#define RECORD_2dL_point record_t<2, double, unsigned long long >::point_t
#define RECORD_2dL_coord record_t<2, double, unsigned long long >::coord_t
#define KDTREE_TYPE_2dL KDTree::KDTree<2, RECORD_2dL, std::pointer_to_binary_function<RECORD_2dL,int,double> >
#define PYKDTREE_2dL PyKDTree<2, double, unsigned long long >




inline bool operator==(RECORD_2dL const& A, RECORD_2dL const& B) {
    return A.data == B.data && A.point[0] == B.point[0] && A.point[1] == B.point[1] ;
}

std::ostream& operator<<(std::ostream& out, RECORD_2dL const& T)
{
    return out << '(' << T.point[0] << ',' << T.point[1] << '|' << T.data << ')';
}



// PyTree_2dL
// PYKDTREE_2dL
// RECORD_2dL
// RECORD_2dL_point
// RECORD_2dL_coord


typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  PYKDTREE_2dL *tree;
} PyTree_2dL;


static PyObject* PyTree_2dL_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyTree_2dL *self;
  
  self = (PyTree_2dL*)type->tp_alloc(type, 0);
  if(self)
    self->tree = NULL;
  
  return (PyObject*)self;
}


static int PyTree_2dL_init(PyTree_2dL *self, PyObject *args, PyObject *kwds)
{
  self->tree = (PYKDTREE_2dL *) new PYKDTREE_2dL();
  return 0;
}


// New datatype dealloc method
static void PyTree_2dL_dealloc(PyTree_2dL* self) {
  if(self->tree)
    delete self->tree;
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PyTree_2dL_repr(PyTree_2dL* self) {

  return PyString_FromString( (char*) "<KDTree_2Double instance>" );
}

static PyObject* PyTree_2dL_str(PyTree_2dL* self) {

  std::ostringstream output;
  std::vector<RECORD_2dL> *temp = 0;
  size_t size;
  size_t full = 6, half = 3;

  if(self && self->tree) {

    size = (int) self->tree->size();

    temp = (std::vector<RECORD_2dL > *) self->tree->get_all();
    
    if (size > full) {
      
      std::vector<RECORD_2dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < half; i++, iter++)
	output << (*iter) << std::endl;
      
      output << "..." << std::endl;
      
      for (size_t i = size-half-1; i < size; i++) 
	output << temp->at(i) << std::endl;
      
    }

    else if( size <= full ) {
    
      std::vector<RECORD_2dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < full; i++, iter++)
        output << (*iter) << std::endl;
    
    }
  
    else {

      output << "[]";
    }
    
    // delete the pointer after use!
    delete temp;
  }

  return PyString_FromString( output.str().c_str() );
}



static PyObject * PyTree_2dL_add(PyTree_2dL* self, PyObject *args, PyObject *kwds) {

  RECORD_2dL record;
  PYKDTREE_2dL *tree;

  if ( !PyArg_ParseTuple(args, "((dd)L)",
			 &record.point[0],
&record.point[1],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (2 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && ((PyTree_2dL*)self)->tree) {

    tree = (PYKDTREE_2dL *)((PyTree_2dL*)self)->tree;
    tree->add(record);

  } else {
    PyErr_SetString(PyExc_TypeError,"Adding record failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject * PyTree_2dL_remove(PyTree_2dL* self, PyObject *args, PyObject *kwds) {

  RECORD_2dL record;
  bool success;
  PyObject *result = Py_False;

  if ( !PyArg_ParseTuple(args, "((dd)L)",
			 &record.point[0],
&record.point[1],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (2 dim double vector, unsigned long long value)");
    return NULL;
  }

  if(self && self->tree) {
    success = self->tree->remove(record);
  } else {
    PyErr_SetString(PyExc_TypeError,"Removing record failed!");
    return NULL;
  }

  if(success)
    result = Py_True;

  Py_INCREF(result);
  return result;  
}

static PyObject * PyTree_2dL_size(PyTree_2dL* self) {

  int result;

  if(self && self->tree) {
    result = (int) self->tree->size();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing size() failed!");
    return NULL;
  }

  return Py_BuildValue((char*)"i", result);
}

static PyObject * PyTree_2dL_optimize(PyTree_2dL* self) {

  if(self && self->tree) {
    self->tree->optimize();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject * PyTree_2dL_find_exact(PyTree_2dL* self, PyObject *args, PyObject *kwds) {

  RECORD_2dL record;
  RECORD_2dL *temp = NULL;
  PyObject *result = NULL;

  if ( !PyArg_ParseTuple(args, "((dd)L)",
			 &record.point[0],
&record.point[1],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (2 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && self->tree) {

    temp = (RECORD_2dL *) self->tree->find_exact(record);

  } else {

    PyErr_SetString(PyExc_TypeError,"find exact failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(dd)",
						 temp->point[0],temp->point[1]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}



static PyObject * PyTree_2dL_find_nearest(PyTree_2dL* self, PyObject *args, PyObject *kwds) {

  RECORD_2dL_point point_t;
  RECORD_2dL_coord *coord_t = NULL;
  RECORD_2dL *temp = NULL;
  PyObject *result = NULL;
  PYKDTREE_2dL *tree;

  if ( !PyArg_ParseTuple(args, "(dd)",
			 &point_t[0],
&point_t[1]
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: 2 dim double vector");
    return NULL;
  }
  
  if(self && ((PyTree_2dL*)self)->tree) {

    tree = (PYKDTREE_2dL *)((PyTree_2dL*)self)->tree;

    coord_t = point_t;

    temp = (RECORD_2dL *) tree->find_nearest(coord_t);

  } else {

    PyErr_SetString(PyExc_TypeError,"find nearest failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(dd)",
						 temp->point[0],temp->point[1]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}

static PyObject * PyTree_2dL_count_within_range(PyTree_2dL* self, PyObject *args, PyObject *kwds) {

  RECORD_2dL_point point_t;
  RECORD_2dL_coord *coord_t = NULL;
  size_t temp;
  PyObject *result = NULL;
  PYKDTREE_2dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(dd)d",
			 &point_t[0],
&point_t[1],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (2 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_2dL*)self)->tree) {

    tree = (PYKDTREE_2dL *)((PyTree_2dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (size_t) tree->count_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"count within range failed!");
    return NULL;
  }

  result = Py_BuildValue("i", (int)temp);
  return result;
}



static PyObject * PyTree_2dL_find_within_range(PyTree_2dL* self, PyObject *args, PyObject *kwds) {

  RECORD_2dL_point point_t;
  RECORD_2dL_coord *coord_t = NULL;
  std::vector<RECORD_2dL> *temp = 0;
  PYKDTREE_2dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(dd)d",
			 &point_t[0],
&point_t[1],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (2 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_2dL*)self)->tree) {

    tree = (PYKDTREE_2dL *)((PyTree_2dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (std::vector<RECORD_2dL > *) tree->find_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"find within range failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_2dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(dd)L",
	(*iter).point[0],(*iter).point[1], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;
  
  return result;
}


static PyObject * PyTree_2dL_get_all(PyTree_2dL* self) {

  std::vector<RECORD_2dL> *temp = 0;

  if(self && self->tree) {

    temp = (std::vector<RECORD_2dL > *) self->tree->get_all();

  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_2dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(dd)L",
	(*iter).point[0],(*iter).point[1], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;

  return result;
}




// Attributes of class are inserted and 'named' here  

static PyMemberDef PyTree_2dL_members[] = {
  // {(char*)"name", T_OBJECT_EX, offsetof(PyTree_2dL, name), 0,
  //  (char*)"name"},

  // {(char*)"datatype", T_INT, offsetof(PyTree_2dL, datatype), 0,
  //  (char*)"datatype"},

  // {(char*)"unit", T_OBJECT_EX, offsetof(PyTree_2dL, unit), 0,
  //  (char*)"unit"},

  // {(char*)"desc", T_OBJECT_EX, offsetof(PyTree_2dL, desc), 0,
  //  (char*)"desc"},
  {NULL}  // Sentinel 
};

// Methods of class are inserted and 'named' here  

// Notice that the function signature type is mapped with a flag 
// METH_NOARGS, METH_VARARGS, or METH_VARARGS | METH_KEYWORDS 

static PyMethodDef PyTree_2dL_methods[] = {

  { (char*) "size",(PyCFunction)PyTree_2dL_size, METH_NOARGS, (char*) ""},

  { (char*) "optimize",(PyCFunction)PyTree_2dL_optimize, METH_NOARGS, (char*) ""},

  { (char*) "get_all",(PyCFunction)PyTree_2dL_get_all, METH_NOARGS, (char*) ""},

  { (char*) "add",(PyCFunction)PyTree_2dL_add, METH_VARARGS, (char*) ""},

  { (char*) "find_exact",(PyCFunction)PyTree_2dL_find_exact, METH_VARARGS, (char*) ""},

  { (char*) "find_nearest",(PyCFunction)PyTree_2dL_find_nearest, METH_VARARGS, (char*) ""},

  { (char*) "find_within_range",(PyCFunction)PyTree_2dL_find_within_range, METH_VARARGS, (char*) ""},

  { (char*) "count_within_range",(PyCFunction)PyTree_2dL_count_within_range, METH_VARARGS, (char*) ""},

  { (char*) "remove",(PyCFunction)PyTree_2dL_remove, METH_VARARGS, (char*) ""},

  {NULL, NULL, 0, NULL}  // Sentinel 
};


static PyTypeObject PyTree_2dLType = {
    PyObject_HEAD_INIT(NULL)
    0,                                // ob_size
    (char*)"KDTree_2Double",           // tp_name, __name__
    sizeof(PyTree_2dL),           // tp_basicsize
    0,                                // tp_itemsize
    (destructor)PyTree_2dL_dealloc, // tp_dealloc
    0,                                // tp_print
    0,                                // tp_getattr, __getattr__
    0,                                // tp_setattr, __setattr__
    0,                                // tp_compare
    (reprfunc)PyTree_2dL_repr,    // tp_repr, __repr__
    0,                                // tp_as_number
    0,                                // tp_as_sequence
    0,                                // tp_as_mapping
    0,                                // tp_hash
    0,                                // tp_call, __call__
    (reprfunc)PyTree_2dL_str,     // tp_str, __str__
    0,                                // tp_getattro
    0,                                // tp_setattro
    0,                                // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags
    (char*)"PyTree_2dL objects",  // tp_doc, __doc__
    0,		                      // tp_traverse 
    0,		                      // tp_clear 
    0,		                      // tp_richcompare 
    0,		                      // tp_weaklistoffset 
    0,		                      // tp_iter, __iter__
    0,		                      // tp_iternext 
    PyTree_2dL_methods,           // tp_methods 
    PyTree_2dL_members,           // tp_members 
    0,                                // tp_getset 
    0,                                // tp_base 
    0,                                // tp_dict, __dict__
    0,                                // tp_descr_get 
    0,                                // tp_descr_set 
    0,                                // tp_dictoffset 
    (initproc)PyTree_2dL_init,    // tp_init, __init__
    0,                                // tp_alloc 
    PyTree_2dL_new,               // tp_new 
};



#define RECORD_3dL record_t<3, double, unsigned long long >
#define RECORD_3dL_point record_t<3, double, unsigned long long >::point_t
#define RECORD_3dL_coord record_t<3, double, unsigned long long >::coord_t
#define KDTREE_TYPE_3dL KDTree::KDTree<3, RECORD_3dL, std::pointer_to_binary_function<RECORD_3dL,int,double> >
#define PYKDTREE_3dL PyKDTree<3, double, unsigned long long >




inline bool operator==(RECORD_3dL const& A, RECORD_3dL const& B) {
    return A.data == B.data && A.point[0] == B.point[0] && A.point[1] == B.point[1] && A.point[2] == B.point[2] ;
}

std::ostream& operator<<(std::ostream& out, RECORD_3dL const& T)
{
    return out << '(' << T.point[0] << ',' << T.point[1] << ',' << T.point[2] << '|' << T.data << ')';
}



// PyTree_3dL
// PYKDTREE_3dL
// RECORD_3dL
// RECORD_3dL_point
// RECORD_3dL_coord


typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  PYKDTREE_3dL *tree;
} PyTree_3dL;


static PyObject* PyTree_3dL_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyTree_3dL *self;
  
  self = (PyTree_3dL*)type->tp_alloc(type, 0);
  if(self)
    self->tree = NULL;
  
  return (PyObject*)self;
}


static int PyTree_3dL_init(PyTree_3dL *self, PyObject *args, PyObject *kwds)
{
  self->tree = (PYKDTREE_3dL *) new PYKDTREE_3dL();
  return 0;
}


// New datatype dealloc method
static void PyTree_3dL_dealloc(PyTree_3dL* self) {
  if(self->tree)
    delete self->tree;
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PyTree_3dL_repr(PyTree_3dL* self) {

  return PyString_FromString( (char*) "<KDTree_3Double instance>" );
}

static PyObject* PyTree_3dL_str(PyTree_3dL* self) {

  std::ostringstream output;
  std::vector<RECORD_3dL> *temp = 0;
  size_t size;
  size_t full = 6, half = 3;

  if(self && self->tree) {

    size = (int) self->tree->size();

    temp = (std::vector<RECORD_3dL > *) self->tree->get_all();
    
    if (size > full) {
      
      std::vector<RECORD_3dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < half; i++, iter++)
	output << (*iter) << std::endl;
      
      output << "..." << std::endl;
      
      for (size_t i = size-half-1; i < size; i++) 
	output << temp->at(i) << std::endl;
      
    }

    else if( size <= full ) {
    
      std::vector<RECORD_3dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < full; i++, iter++)
        output << (*iter) << std::endl;
    
    }
  
    else {

      output << "[]";
    }
    
    // delete the pointer after use!
    delete temp;
  }

  return PyString_FromString( output.str().c_str() );
}



static PyObject * PyTree_3dL_add(PyTree_3dL* self, PyObject *args, PyObject *kwds) {

  RECORD_3dL record;
  PYKDTREE_3dL *tree;

  if ( !PyArg_ParseTuple(args, "((ddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (3 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && ((PyTree_3dL*)self)->tree) {

    tree = (PYKDTREE_3dL *)((PyTree_3dL*)self)->tree;
    tree->add(record);

  } else {
    PyErr_SetString(PyExc_TypeError,"Adding record failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject * PyTree_3dL_remove(PyTree_3dL* self, PyObject *args, PyObject *kwds) {

  RECORD_3dL record;
  bool success;
  PyObject *result = Py_False;

  if ( !PyArg_ParseTuple(args, "((ddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (3 dim double vector, unsigned long long value)");
    return NULL;
  }

  if(self && self->tree) {
    success = self->tree->remove(record);
  } else {
    PyErr_SetString(PyExc_TypeError,"Removing record failed!");
    return NULL;
  }

  if(success)
    result = Py_True;

  Py_INCREF(result);
  return result;  
}

static PyObject * PyTree_3dL_size(PyTree_3dL* self) {

  int result;

  if(self && self->tree) {
    result = (int) self->tree->size();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing size() failed!");
    return NULL;
  }

  return Py_BuildValue((char*)"i", result);
}

static PyObject * PyTree_3dL_optimize(PyTree_3dL* self) {

  if(self && self->tree) {
    self->tree->optimize();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject * PyTree_3dL_find_exact(PyTree_3dL* self, PyObject *args, PyObject *kwds) {

  RECORD_3dL record;
  RECORD_3dL *temp = NULL;
  PyObject *result = NULL;

  if ( !PyArg_ParseTuple(args, "((ddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (3 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && self->tree) {

    temp = (RECORD_3dL *) self->tree->find_exact(record);

  } else {

    PyErr_SetString(PyExc_TypeError,"find exact failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(ddd)",
						 temp->point[0],temp->point[1],temp->point[2]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}



static PyObject * PyTree_3dL_find_nearest(PyTree_3dL* self, PyObject *args, PyObject *kwds) {

  RECORD_3dL_point point_t;
  RECORD_3dL_coord *coord_t = NULL;
  RECORD_3dL *temp = NULL;
  PyObject *result = NULL;
  PYKDTREE_3dL *tree;

  if ( !PyArg_ParseTuple(args, "(ddd)",
			 &point_t[0],
&point_t[1],
&point_t[2]
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 3 elements: 3 dim double vector");
    return NULL;
  }
  
  if(self && ((PyTree_3dL*)self)->tree) {

    tree = (PYKDTREE_3dL *)((PyTree_3dL*)self)->tree;

    coord_t = point_t;

    temp = (RECORD_3dL *) tree->find_nearest(coord_t);

  } else {

    PyErr_SetString(PyExc_TypeError,"find nearest failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(ddd)",
						 temp->point[0],temp->point[1],temp->point[2]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}

static PyObject * PyTree_3dL_count_within_range(PyTree_3dL* self, PyObject *args, PyObject *kwds) {

  RECORD_3dL_point point_t;
  RECORD_3dL_coord *coord_t = NULL;
  size_t temp;
  PyObject *result = NULL;
  PYKDTREE_3dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(ddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (3 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_3dL*)self)->tree) {

    tree = (PYKDTREE_3dL *)((PyTree_3dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (size_t) tree->count_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"count within range failed!");
    return NULL;
  }

  result = Py_BuildValue("i", (int)temp);
  return result;
}



static PyObject * PyTree_3dL_find_within_range(PyTree_3dL* self, PyObject *args, PyObject *kwds) {

  RECORD_3dL_point point_t;
  RECORD_3dL_coord *coord_t = NULL;
  std::vector<RECORD_3dL> *temp = 0;
  PYKDTREE_3dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(ddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (3 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_3dL*)self)->tree) {

    tree = (PYKDTREE_3dL *)((PyTree_3dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (std::vector<RECORD_3dL > *) tree->find_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"find within range failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_3dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(ddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;
  
  return result;
}


static PyObject * PyTree_3dL_get_all(PyTree_3dL* self) {

  std::vector<RECORD_3dL> *temp = 0;

  if(self && self->tree) {

    temp = (std::vector<RECORD_3dL > *) self->tree->get_all();

  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_3dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(ddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;

  return result;
}




// Attributes of class are inserted and 'named' here  

static PyMemberDef PyTree_3dL_members[] = {
  // {(char*)"name", T_OBJECT_EX, offsetof(PyTree_3dL, name), 0,
  //  (char*)"name"},

  // {(char*)"datatype", T_INT, offsetof(PyTree_3dL, datatype), 0,
  //  (char*)"datatype"},

  // {(char*)"unit", T_OBJECT_EX, offsetof(PyTree_3dL, unit), 0,
  //  (char*)"unit"},

  // {(char*)"desc", T_OBJECT_EX, offsetof(PyTree_3dL, desc), 0,
  //  (char*)"desc"},
  {NULL}  // Sentinel 
};

// Methods of class are inserted and 'named' here  

// Notice that the function signature type is mapped with a flag 
// METH_NOARGS, METH_VARARGS, or METH_VARARGS | METH_KEYWORDS 

static PyMethodDef PyTree_3dL_methods[] = {

  { (char*) "size",(PyCFunction)PyTree_3dL_size, METH_NOARGS, (char*) ""},

  { (char*) "optimize",(PyCFunction)PyTree_3dL_optimize, METH_NOARGS, (char*) ""},

  { (char*) "get_all",(PyCFunction)PyTree_3dL_get_all, METH_NOARGS, (char*) ""},

  { (char*) "add",(PyCFunction)PyTree_3dL_add, METH_VARARGS, (char*) ""},

  { (char*) "find_exact",(PyCFunction)PyTree_3dL_find_exact, METH_VARARGS, (char*) ""},

  { (char*) "find_nearest",(PyCFunction)PyTree_3dL_find_nearest, METH_VARARGS, (char*) ""},

  { (char*) "find_within_range",(PyCFunction)PyTree_3dL_find_within_range, METH_VARARGS, (char*) ""},

  { (char*) "count_within_range",(PyCFunction)PyTree_3dL_count_within_range, METH_VARARGS, (char*) ""},

  { (char*) "remove",(PyCFunction)PyTree_3dL_remove, METH_VARARGS, (char*) ""},

  {NULL, NULL, 0, NULL}  // Sentinel 
};


static PyTypeObject PyTree_3dLType = {
    PyObject_HEAD_INIT(NULL)
    0,                                // ob_size
    (char*)"KDTree_3Double",           // tp_name, __name__
    sizeof(PyTree_3dL),           // tp_basicsize
    0,                                // tp_itemsize
    (destructor)PyTree_3dL_dealloc, // tp_dealloc
    0,                                // tp_print
    0,                                // tp_getattr, __getattr__
    0,                                // tp_setattr, __setattr__
    0,                                // tp_compare
    (reprfunc)PyTree_3dL_repr,    // tp_repr, __repr__
    0,                                // tp_as_number
    0,                                // tp_as_sequence
    0,                                // tp_as_mapping
    0,                                // tp_hash
    0,                                // tp_call, __call__
    (reprfunc)PyTree_3dL_str,     // tp_str, __str__
    0,                                // tp_getattro
    0,                                // tp_setattro
    0,                                // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags
    (char*)"PyTree_3dL objects",  // tp_doc, __doc__
    0,		                      // tp_traverse 
    0,		                      // tp_clear 
    0,		                      // tp_richcompare 
    0,		                      // tp_weaklistoffset 
    0,		                      // tp_iter, __iter__
    0,		                      // tp_iternext 
    PyTree_3dL_methods,           // tp_methods 
    PyTree_3dL_members,           // tp_members 
    0,                                // tp_getset 
    0,                                // tp_base 
    0,                                // tp_dict, __dict__
    0,                                // tp_descr_get 
    0,                                // tp_descr_set 
    0,                                // tp_dictoffset 
    (initproc)PyTree_3dL_init,    // tp_init, __init__
    0,                                // tp_alloc 
    PyTree_3dL_new,               // tp_new 
};



#define RECORD_4dL record_t<4, double, unsigned long long >
#define RECORD_4dL_point record_t<4, double, unsigned long long >::point_t
#define RECORD_4dL_coord record_t<4, double, unsigned long long >::coord_t
#define KDTREE_TYPE_4dL KDTree::KDTree<4, RECORD_4dL, std::pointer_to_binary_function<RECORD_4dL,int,double> >
#define PYKDTREE_4dL PyKDTree<4, double, unsigned long long >




inline bool operator==(RECORD_4dL const& A, RECORD_4dL const& B) {
    return A.data == B.data && A.point[0] == B.point[0] && A.point[1] == B.point[1] && A.point[2] == B.point[2] && A.point[3] == B.point[3] ;
}

std::ostream& operator<<(std::ostream& out, RECORD_4dL const& T)
{
    return out << '(' << T.point[0] << ',' << T.point[1] << ',' << T.point[2] << ',' << T.point[3] << '|' << T.data << ')';
}



// PyTree_4dL
// PYKDTREE_4dL
// RECORD_4dL
// RECORD_4dL_point
// RECORD_4dL_coord


typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  PYKDTREE_4dL *tree;
} PyTree_4dL;


static PyObject* PyTree_4dL_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyTree_4dL *self;
  
  self = (PyTree_4dL*)type->tp_alloc(type, 0);
  if(self)
    self->tree = NULL;
  
  return (PyObject*)self;
}


static int PyTree_4dL_init(PyTree_4dL *self, PyObject *args, PyObject *kwds)
{
  self->tree = (PYKDTREE_4dL *) new PYKDTREE_4dL();
  return 0;
}


// New datatype dealloc method
static void PyTree_4dL_dealloc(PyTree_4dL* self) {
  if(self->tree)
    delete self->tree;
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PyTree_4dL_repr(PyTree_4dL* self) {

  return PyString_FromString( (char*) "<KDTree_4Double instance>" );
}

static PyObject* PyTree_4dL_str(PyTree_4dL* self) {

  std::ostringstream output;
  std::vector<RECORD_4dL> *temp = 0;
  size_t size;
  size_t full = 6, half = 3;

  if(self && self->tree) {

    size = (int) self->tree->size();

    temp = (std::vector<RECORD_4dL > *) self->tree->get_all();
    
    if (size > full) {
      
      std::vector<RECORD_4dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < half; i++, iter++)
	output << (*iter) << std::endl;
      
      output << "..." << std::endl;
      
      for (size_t i = size-half-1; i < size; i++) 
	output << temp->at(i) << std::endl;
      
    }

    else if( size <= full ) {
    
      std::vector<RECORD_4dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < full; i++, iter++)
        output << (*iter) << std::endl;
    
    }
  
    else {

      output << "[]";
    }
    
    // delete the pointer after use!
    delete temp;
  }

  return PyString_FromString( output.str().c_str() );
}



static PyObject * PyTree_4dL_add(PyTree_4dL* self, PyObject *args, PyObject *kwds) {

  RECORD_4dL record;
  PYKDTREE_4dL *tree;

  if ( !PyArg_ParseTuple(args, "((dddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (4 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && ((PyTree_4dL*)self)->tree) {

    tree = (PYKDTREE_4dL *)((PyTree_4dL*)self)->tree;
    tree->add(record);

  } else {
    PyErr_SetString(PyExc_TypeError,"Adding record failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject * PyTree_4dL_remove(PyTree_4dL* self, PyObject *args, PyObject *kwds) {

  RECORD_4dL record;
  bool success;
  PyObject *result = Py_False;

  if ( !PyArg_ParseTuple(args, "((dddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (4 dim double vector, unsigned long long value)");
    return NULL;
  }

  if(self && self->tree) {
    success = self->tree->remove(record);
  } else {
    PyErr_SetString(PyExc_TypeError,"Removing record failed!");
    return NULL;
  }

  if(success)
    result = Py_True;

  Py_INCREF(result);
  return result;  
}

static PyObject * PyTree_4dL_size(PyTree_4dL* self) {

  int result;

  if(self && self->tree) {
    result = (int) self->tree->size();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing size() failed!");
    return NULL;
  }

  return Py_BuildValue((char*)"i", result);
}

static PyObject * PyTree_4dL_optimize(PyTree_4dL* self) {

  if(self && self->tree) {
    self->tree->optimize();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject * PyTree_4dL_find_exact(PyTree_4dL* self, PyObject *args, PyObject *kwds) {

  RECORD_4dL record;
  RECORD_4dL *temp = NULL;
  PyObject *result = NULL;

  if ( !PyArg_ParseTuple(args, "((dddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (4 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && self->tree) {

    temp = (RECORD_4dL *) self->tree->find_exact(record);

  } else {

    PyErr_SetString(PyExc_TypeError,"find exact failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(dddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}



static PyObject * PyTree_4dL_find_nearest(PyTree_4dL* self, PyObject *args, PyObject *kwds) {

  RECORD_4dL_point point_t;
  RECORD_4dL_coord *coord_t = NULL;
  RECORD_4dL *temp = NULL;
  PyObject *result = NULL;
  PYKDTREE_4dL *tree;

  if ( !PyArg_ParseTuple(args, "(dddd)",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3]
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 4 elements: 4 dim double vector");
    return NULL;
  }
  
  if(self && ((PyTree_4dL*)self)->tree) {

    tree = (PYKDTREE_4dL *)((PyTree_4dL*)self)->tree;

    coord_t = point_t;

    temp = (RECORD_4dL *) tree->find_nearest(coord_t);

  } else {

    PyErr_SetString(PyExc_TypeError,"find nearest failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(dddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}

static PyObject * PyTree_4dL_count_within_range(PyTree_4dL* self, PyObject *args, PyObject *kwds) {

  RECORD_4dL_point point_t;
  RECORD_4dL_coord *coord_t = NULL;
  size_t temp;
  PyObject *result = NULL;
  PYKDTREE_4dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(dddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (4 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_4dL*)self)->tree) {

    tree = (PYKDTREE_4dL *)((PyTree_4dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (size_t) tree->count_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"count within range failed!");
    return NULL;
  }

  result = Py_BuildValue("i", (int)temp);
  return result;
}



static PyObject * PyTree_4dL_find_within_range(PyTree_4dL* self, PyObject *args, PyObject *kwds) {

  RECORD_4dL_point point_t;
  RECORD_4dL_coord *coord_t = NULL;
  std::vector<RECORD_4dL> *temp = 0;
  PYKDTREE_4dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(dddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (4 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_4dL*)self)->tree) {

    tree = (PYKDTREE_4dL *)((PyTree_4dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (std::vector<RECORD_4dL > *) tree->find_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"find within range failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_4dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(dddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;
  
  return result;
}


static PyObject * PyTree_4dL_get_all(PyTree_4dL* self) {

  std::vector<RECORD_4dL> *temp = 0;

  if(self && self->tree) {

    temp = (std::vector<RECORD_4dL > *) self->tree->get_all();

  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_4dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(dddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;

  return result;
}




// Attributes of class are inserted and 'named' here  

static PyMemberDef PyTree_4dL_members[] = {
  // {(char*)"name", T_OBJECT_EX, offsetof(PyTree_4dL, name), 0,
  //  (char*)"name"},

  // {(char*)"datatype", T_INT, offsetof(PyTree_4dL, datatype), 0,
  //  (char*)"datatype"},

  // {(char*)"unit", T_OBJECT_EX, offsetof(PyTree_4dL, unit), 0,
  //  (char*)"unit"},

  // {(char*)"desc", T_OBJECT_EX, offsetof(PyTree_4dL, desc), 0,
  //  (char*)"desc"},
  {NULL}  // Sentinel 
};

// Methods of class are inserted and 'named' here  

// Notice that the function signature type is mapped with a flag 
// METH_NOARGS, METH_VARARGS, or METH_VARARGS | METH_KEYWORDS 

static PyMethodDef PyTree_4dL_methods[] = {

  { (char*) "size",(PyCFunction)PyTree_4dL_size, METH_NOARGS, (char*) ""},

  { (char*) "optimize",(PyCFunction)PyTree_4dL_optimize, METH_NOARGS, (char*) ""},

  { (char*) "get_all",(PyCFunction)PyTree_4dL_get_all, METH_NOARGS, (char*) ""},

  { (char*) "add",(PyCFunction)PyTree_4dL_add, METH_VARARGS, (char*) ""},

  { (char*) "find_exact",(PyCFunction)PyTree_4dL_find_exact, METH_VARARGS, (char*) ""},

  { (char*) "find_nearest",(PyCFunction)PyTree_4dL_find_nearest, METH_VARARGS, (char*) ""},

  { (char*) "find_within_range",(PyCFunction)PyTree_4dL_find_within_range, METH_VARARGS, (char*) ""},

  { (char*) "count_within_range",(PyCFunction)PyTree_4dL_count_within_range, METH_VARARGS, (char*) ""},

  { (char*) "remove",(PyCFunction)PyTree_4dL_remove, METH_VARARGS, (char*) ""},

  {NULL, NULL, 0, NULL}  // Sentinel 
};


static PyTypeObject PyTree_4dLType = {
    PyObject_HEAD_INIT(NULL)
    0,                                // ob_size
    (char*)"KDTree_4Double",           // tp_name, __name__
    sizeof(PyTree_4dL),           // tp_basicsize
    0,                                // tp_itemsize
    (destructor)PyTree_4dL_dealloc, // tp_dealloc
    0,                                // tp_print
    0,                                // tp_getattr, __getattr__
    0,                                // tp_setattr, __setattr__
    0,                                // tp_compare
    (reprfunc)PyTree_4dL_repr,    // tp_repr, __repr__
    0,                                // tp_as_number
    0,                                // tp_as_sequence
    0,                                // tp_as_mapping
    0,                                // tp_hash
    0,                                // tp_call, __call__
    (reprfunc)PyTree_4dL_str,     // tp_str, __str__
    0,                                // tp_getattro
    0,                                // tp_setattro
    0,                                // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags
    (char*)"PyTree_4dL objects",  // tp_doc, __doc__
    0,		                      // tp_traverse 
    0,		                      // tp_clear 
    0,		                      // tp_richcompare 
    0,		                      // tp_weaklistoffset 
    0,		                      // tp_iter, __iter__
    0,		                      // tp_iternext 
    PyTree_4dL_methods,           // tp_methods 
    PyTree_4dL_members,           // tp_members 
    0,                                // tp_getset 
    0,                                // tp_base 
    0,                                // tp_dict, __dict__
    0,                                // tp_descr_get 
    0,                                // tp_descr_set 
    0,                                // tp_dictoffset 
    (initproc)PyTree_4dL_init,    // tp_init, __init__
    0,                                // tp_alloc 
    PyTree_4dL_new,               // tp_new 
};



#define RECORD_5dL record_t<5, double, unsigned long long >
#define RECORD_5dL_point record_t<5, double, unsigned long long >::point_t
#define RECORD_5dL_coord record_t<5, double, unsigned long long >::coord_t
#define KDTREE_TYPE_5dL KDTree::KDTree<5, RECORD_5dL, std::pointer_to_binary_function<RECORD_5dL,int,double> >
#define PYKDTREE_5dL PyKDTree<5, double, unsigned long long >




inline bool operator==(RECORD_5dL const& A, RECORD_5dL const& B) {
    return A.data == B.data && A.point[0] == B.point[0] && A.point[1] == B.point[1] && A.point[2] == B.point[2] && A.point[3] == B.point[3] && A.point[4] == B.point[4] ;
}

std::ostream& operator<<(std::ostream& out, RECORD_5dL const& T)
{
    return out << '(' << T.point[0] << ',' << T.point[1] << ',' << T.point[2] << ',' << T.point[3] << ',' << T.point[4] << '|' << T.data << ')';
}



// PyTree_5dL
// PYKDTREE_5dL
// RECORD_5dL
// RECORD_5dL_point
// RECORD_5dL_coord


typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  PYKDTREE_5dL *tree;
} PyTree_5dL;


static PyObject* PyTree_5dL_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyTree_5dL *self;
  
  self = (PyTree_5dL*)type->tp_alloc(type, 0);
  if(self)
    self->tree = NULL;
  
  return (PyObject*)self;
}


static int PyTree_5dL_init(PyTree_5dL *self, PyObject *args, PyObject *kwds)
{
  self->tree = (PYKDTREE_5dL *) new PYKDTREE_5dL();
  return 0;
}


// New datatype dealloc method
static void PyTree_5dL_dealloc(PyTree_5dL* self) {
  if(self->tree)
    delete self->tree;
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PyTree_5dL_repr(PyTree_5dL* self) {

  return PyString_FromString( (char*) "<KDTree_5Double instance>" );
}

static PyObject* PyTree_5dL_str(PyTree_5dL* self) {

  std::ostringstream output;
  std::vector<RECORD_5dL> *temp = 0;
  size_t size;
  size_t full = 6, half = 3;

  if(self && self->tree) {

    size = (int) self->tree->size();

    temp = (std::vector<RECORD_5dL > *) self->tree->get_all();
    
    if (size > full) {
      
      std::vector<RECORD_5dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < half; i++, iter++)
	output << (*iter) << std::endl;
      
      output << "..." << std::endl;
      
      for (size_t i = size-half-1; i < size; i++) 
	output << temp->at(i) << std::endl;
      
    }

    else if( size <= full ) {
    
      std::vector<RECORD_5dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < full; i++, iter++)
        output << (*iter) << std::endl;
    
    }
  
    else {

      output << "[]";
    }
    
    // delete the pointer after use!
    delete temp;
  }

  return PyString_FromString( output.str().c_str() );
}



static PyObject * PyTree_5dL_add(PyTree_5dL* self, PyObject *args, PyObject *kwds) {

  RECORD_5dL record;
  PYKDTREE_5dL *tree;

  if ( !PyArg_ParseTuple(args, "((ddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (5 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && ((PyTree_5dL*)self)->tree) {

    tree = (PYKDTREE_5dL *)((PyTree_5dL*)self)->tree;
    tree->add(record);

  } else {
    PyErr_SetString(PyExc_TypeError,"Adding record failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject * PyTree_5dL_remove(PyTree_5dL* self, PyObject *args, PyObject *kwds) {

  RECORD_5dL record;
  bool success;
  PyObject *result = Py_False;

  if ( !PyArg_ParseTuple(args, "((ddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (5 dim double vector, unsigned long long value)");
    return NULL;
  }

  if(self && self->tree) {
    success = self->tree->remove(record);
  } else {
    PyErr_SetString(PyExc_TypeError,"Removing record failed!");
    return NULL;
  }

  if(success)
    result = Py_True;

  Py_INCREF(result);
  return result;  
}

static PyObject * PyTree_5dL_size(PyTree_5dL* self) {

  int result;

  if(self && self->tree) {
    result = (int) self->tree->size();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing size() failed!");
    return NULL;
  }

  return Py_BuildValue((char*)"i", result);
}

static PyObject * PyTree_5dL_optimize(PyTree_5dL* self) {

  if(self && self->tree) {
    self->tree->optimize();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject * PyTree_5dL_find_exact(PyTree_5dL* self, PyObject *args, PyObject *kwds) {

  RECORD_5dL record;
  RECORD_5dL *temp = NULL;
  PyObject *result = NULL;

  if ( !PyArg_ParseTuple(args, "((ddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (5 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && self->tree) {

    temp = (RECORD_5dL *) self->tree->find_exact(record);

  } else {

    PyErr_SetString(PyExc_TypeError,"find exact failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(ddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}



static PyObject * PyTree_5dL_find_nearest(PyTree_5dL* self, PyObject *args, PyObject *kwds) {

  RECORD_5dL_point point_t;
  RECORD_5dL_coord *coord_t = NULL;
  RECORD_5dL *temp = NULL;
  PyObject *result = NULL;
  PYKDTREE_5dL *tree;

  if ( !PyArg_ParseTuple(args, "(ddddd)",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4]
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 5 elements: 5 dim double vector");
    return NULL;
  }
  
  if(self && ((PyTree_5dL*)self)->tree) {

    tree = (PYKDTREE_5dL *)((PyTree_5dL*)self)->tree;

    coord_t = point_t;

    temp = (RECORD_5dL *) tree->find_nearest(coord_t);

  } else {

    PyErr_SetString(PyExc_TypeError,"find nearest failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(ddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}

static PyObject * PyTree_5dL_count_within_range(PyTree_5dL* self, PyObject *args, PyObject *kwds) {

  RECORD_5dL_point point_t;
  RECORD_5dL_coord *coord_t = NULL;
  size_t temp;
  PyObject *result = NULL;
  PYKDTREE_5dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(ddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (5 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_5dL*)self)->tree) {

    tree = (PYKDTREE_5dL *)((PyTree_5dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (size_t) tree->count_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"count within range failed!");
    return NULL;
  }

  result = Py_BuildValue("i", (int)temp);
  return result;
}



static PyObject * PyTree_5dL_find_within_range(PyTree_5dL* self, PyObject *args, PyObject *kwds) {

  RECORD_5dL_point point_t;
  RECORD_5dL_coord *coord_t = NULL;
  std::vector<RECORD_5dL> *temp = 0;
  PYKDTREE_5dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(ddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (5 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_5dL*)self)->tree) {

    tree = (PYKDTREE_5dL *)((PyTree_5dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (std::vector<RECORD_5dL > *) tree->find_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"find within range failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_5dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(ddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;
  
  return result;
}


static PyObject * PyTree_5dL_get_all(PyTree_5dL* self) {

  std::vector<RECORD_5dL> *temp = 0;

  if(self && self->tree) {

    temp = (std::vector<RECORD_5dL > *) self->tree->get_all();

  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_5dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(ddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;

  return result;
}




// Attributes of class are inserted and 'named' here  

static PyMemberDef PyTree_5dL_members[] = {
  // {(char*)"name", T_OBJECT_EX, offsetof(PyTree_5dL, name), 0,
  //  (char*)"name"},

  // {(char*)"datatype", T_INT, offsetof(PyTree_5dL, datatype), 0,
  //  (char*)"datatype"},

  // {(char*)"unit", T_OBJECT_EX, offsetof(PyTree_5dL, unit), 0,
  //  (char*)"unit"},

  // {(char*)"desc", T_OBJECT_EX, offsetof(PyTree_5dL, desc), 0,
  //  (char*)"desc"},
  {NULL}  // Sentinel 
};

// Methods of class are inserted and 'named' here  

// Notice that the function signature type is mapped with a flag 
// METH_NOARGS, METH_VARARGS, or METH_VARARGS | METH_KEYWORDS 

static PyMethodDef PyTree_5dL_methods[] = {

  { (char*) "size",(PyCFunction)PyTree_5dL_size, METH_NOARGS, (char*) ""},

  { (char*) "optimize",(PyCFunction)PyTree_5dL_optimize, METH_NOARGS, (char*) ""},

  { (char*) "get_all",(PyCFunction)PyTree_5dL_get_all, METH_NOARGS, (char*) ""},

  { (char*) "add",(PyCFunction)PyTree_5dL_add, METH_VARARGS, (char*) ""},

  { (char*) "find_exact",(PyCFunction)PyTree_5dL_find_exact, METH_VARARGS, (char*) ""},

  { (char*) "find_nearest",(PyCFunction)PyTree_5dL_find_nearest, METH_VARARGS, (char*) ""},

  { (char*) "find_within_range",(PyCFunction)PyTree_5dL_find_within_range, METH_VARARGS, (char*) ""},

  { (char*) "count_within_range",(PyCFunction)PyTree_5dL_count_within_range, METH_VARARGS, (char*) ""},

  { (char*) "remove",(PyCFunction)PyTree_5dL_remove, METH_VARARGS, (char*) ""},

  {NULL, NULL, 0, NULL}  // Sentinel 
};


static PyTypeObject PyTree_5dLType = {
    PyObject_HEAD_INIT(NULL)
    0,                                // ob_size
    (char*)"KDTree_5Double",           // tp_name, __name__
    sizeof(PyTree_5dL),           // tp_basicsize
    0,                                // tp_itemsize
    (destructor)PyTree_5dL_dealloc, // tp_dealloc
    0,                                // tp_print
    0,                                // tp_getattr, __getattr__
    0,                                // tp_setattr, __setattr__
    0,                                // tp_compare
    (reprfunc)PyTree_5dL_repr,    // tp_repr, __repr__
    0,                                // tp_as_number
    0,                                // tp_as_sequence
    0,                                // tp_as_mapping
    0,                                // tp_hash
    0,                                // tp_call, __call__
    (reprfunc)PyTree_5dL_str,     // tp_str, __str__
    0,                                // tp_getattro
    0,                                // tp_setattro
    0,                                // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags
    (char*)"PyTree_5dL objects",  // tp_doc, __doc__
    0,		                      // tp_traverse 
    0,		                      // tp_clear 
    0,		                      // tp_richcompare 
    0,		                      // tp_weaklistoffset 
    0,		                      // tp_iter, __iter__
    0,		                      // tp_iternext 
    PyTree_5dL_methods,           // tp_methods 
    PyTree_5dL_members,           // tp_members 
    0,                                // tp_getset 
    0,                                // tp_base 
    0,                                // tp_dict, __dict__
    0,                                // tp_descr_get 
    0,                                // tp_descr_set 
    0,                                // tp_dictoffset 
    (initproc)PyTree_5dL_init,    // tp_init, __init__
    0,                                // tp_alloc 
    PyTree_5dL_new,               // tp_new 
};



#define RECORD_6dL record_t<6, double, unsigned long long >
#define RECORD_6dL_point record_t<6, double, unsigned long long >::point_t
#define RECORD_6dL_coord record_t<6, double, unsigned long long >::coord_t
#define KDTREE_TYPE_6dL KDTree::KDTree<6, RECORD_6dL, std::pointer_to_binary_function<RECORD_6dL,int,double> >
#define PYKDTREE_6dL PyKDTree<6, double, unsigned long long >




inline bool operator==(RECORD_6dL const& A, RECORD_6dL const& B) {
    return A.data == B.data && A.point[0] == B.point[0] && A.point[1] == B.point[1] && A.point[2] == B.point[2] && A.point[3] == B.point[3] && A.point[4] == B.point[4] && A.point[5] == B.point[5] ;
}

std::ostream& operator<<(std::ostream& out, RECORD_6dL const& T)
{
    return out << '(' << T.point[0] << ',' << T.point[1] << ',' << T.point[2] << ',' << T.point[3] << ',' << T.point[4] << ',' << T.point[5] << '|' << T.data << ')';
}



// PyTree_6dL
// PYKDTREE_6dL
// RECORD_6dL
// RECORD_6dL_point
// RECORD_6dL_coord


typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  PYKDTREE_6dL *tree;
} PyTree_6dL;


static PyObject* PyTree_6dL_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyTree_6dL *self;
  
  self = (PyTree_6dL*)type->tp_alloc(type, 0);
  if(self)
    self->tree = NULL;
  
  return (PyObject*)self;
}


static int PyTree_6dL_init(PyTree_6dL *self, PyObject *args, PyObject *kwds)
{
  self->tree = (PYKDTREE_6dL *) new PYKDTREE_6dL();
  return 0;
}


// New datatype dealloc method
static void PyTree_6dL_dealloc(PyTree_6dL* self) {
  if(self->tree)
    delete self->tree;
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PyTree_6dL_repr(PyTree_6dL* self) {

  return PyString_FromString( (char*) "<KDTree_6Double instance>" );
}

static PyObject* PyTree_6dL_str(PyTree_6dL* self) {

  std::ostringstream output;
  std::vector<RECORD_6dL> *temp = 0;
  size_t size;
  size_t full = 6, half = 3;

  if(self && self->tree) {

    size = (int) self->tree->size();

    temp = (std::vector<RECORD_6dL > *) self->tree->get_all();
    
    if (size > full) {
      
      std::vector<RECORD_6dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < half; i++, iter++)
	output << (*iter) << std::endl;
      
      output << "..." << std::endl;
      
      for (size_t i = size-half-1; i < size; i++) 
	output << temp->at(i) << std::endl;
      
    }

    else if( size <= full ) {
    
      std::vector<RECORD_6dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < full; i++, iter++)
        output << (*iter) << std::endl;
    
    }
  
    else {

      output << "[]";
    }
    
    // delete the pointer after use!
    delete temp;
  }

  return PyString_FromString( output.str().c_str() );
}



static PyObject * PyTree_6dL_add(PyTree_6dL* self, PyObject *args, PyObject *kwds) {

  RECORD_6dL record;
  PYKDTREE_6dL *tree;

  if ( !PyArg_ParseTuple(args, "((dddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (6 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && ((PyTree_6dL*)self)->tree) {

    tree = (PYKDTREE_6dL *)((PyTree_6dL*)self)->tree;
    tree->add(record);

  } else {
    PyErr_SetString(PyExc_TypeError,"Adding record failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject * PyTree_6dL_remove(PyTree_6dL* self, PyObject *args, PyObject *kwds) {

  RECORD_6dL record;
  bool success;
  PyObject *result = Py_False;

  if ( !PyArg_ParseTuple(args, "((dddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (6 dim double vector, unsigned long long value)");
    return NULL;
  }

  if(self && self->tree) {
    success = self->tree->remove(record);
  } else {
    PyErr_SetString(PyExc_TypeError,"Removing record failed!");
    return NULL;
  }

  if(success)
    result = Py_True;

  Py_INCREF(result);
  return result;  
}

static PyObject * PyTree_6dL_size(PyTree_6dL* self) {

  int result;

  if(self && self->tree) {
    result = (int) self->tree->size();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing size() failed!");
    return NULL;
  }

  return Py_BuildValue((char*)"i", result);
}

static PyObject * PyTree_6dL_optimize(PyTree_6dL* self) {

  if(self && self->tree) {
    self->tree->optimize();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject * PyTree_6dL_find_exact(PyTree_6dL* self, PyObject *args, PyObject *kwds) {

  RECORD_6dL record;
  RECORD_6dL *temp = NULL;
  PyObject *result = NULL;

  if ( !PyArg_ParseTuple(args, "((dddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (6 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && self->tree) {

    temp = (RECORD_6dL *) self->tree->find_exact(record);

  } else {

    PyErr_SetString(PyExc_TypeError,"find exact failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(dddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4],temp->point[5]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}



static PyObject * PyTree_6dL_find_nearest(PyTree_6dL* self, PyObject *args, PyObject *kwds) {

  RECORD_6dL_point point_t;
  RECORD_6dL_coord *coord_t = NULL;
  RECORD_6dL *temp = NULL;
  PyObject *result = NULL;
  PYKDTREE_6dL *tree;

  if ( !PyArg_ParseTuple(args, "(dddddd)",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5]
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 6 elements: 6 dim double vector");
    return NULL;
  }
  
  if(self && ((PyTree_6dL*)self)->tree) {

    tree = (PYKDTREE_6dL *)((PyTree_6dL*)self)->tree;

    coord_t = point_t;

    temp = (RECORD_6dL *) tree->find_nearest(coord_t);

  } else {

    PyErr_SetString(PyExc_TypeError,"find nearest failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(dddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4],temp->point[5]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}

static PyObject * PyTree_6dL_count_within_range(PyTree_6dL* self, PyObject *args, PyObject *kwds) {

  RECORD_6dL_point point_t;
  RECORD_6dL_coord *coord_t = NULL;
  size_t temp;
  PyObject *result = NULL;
  PYKDTREE_6dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(dddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (6 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_6dL*)self)->tree) {

    tree = (PYKDTREE_6dL *)((PyTree_6dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (size_t) tree->count_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"count within range failed!");
    return NULL;
  }

  result = Py_BuildValue("i", (int)temp);
  return result;
}



static PyObject * PyTree_6dL_find_within_range(PyTree_6dL* self, PyObject *args, PyObject *kwds) {

  RECORD_6dL_point point_t;
  RECORD_6dL_coord *coord_t = NULL;
  std::vector<RECORD_6dL> *temp = 0;
  PYKDTREE_6dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(dddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (6 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_6dL*)self)->tree) {

    tree = (PYKDTREE_6dL *)((PyTree_6dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (std::vector<RECORD_6dL > *) tree->find_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"find within range failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_6dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(dddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4],(*iter).point[5], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;
  
  return result;
}


static PyObject * PyTree_6dL_get_all(PyTree_6dL* self) {

  std::vector<RECORD_6dL> *temp = 0;

  if(self && self->tree) {

    temp = (std::vector<RECORD_6dL > *) self->tree->get_all();

  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_6dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(dddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4],(*iter).point[5], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;

  return result;
}




// Attributes of class are inserted and 'named' here  

static PyMemberDef PyTree_6dL_members[] = {
  // {(char*)"name", T_OBJECT_EX, offsetof(PyTree_6dL, name), 0,
  //  (char*)"name"},

  // {(char*)"datatype", T_INT, offsetof(PyTree_6dL, datatype), 0,
  //  (char*)"datatype"},

  // {(char*)"unit", T_OBJECT_EX, offsetof(PyTree_6dL, unit), 0,
  //  (char*)"unit"},

  // {(char*)"desc", T_OBJECT_EX, offsetof(PyTree_6dL, desc), 0,
  //  (char*)"desc"},
  {NULL}  // Sentinel 
};

// Methods of class are inserted and 'named' here  

// Notice that the function signature type is mapped with a flag 
// METH_NOARGS, METH_VARARGS, or METH_VARARGS | METH_KEYWORDS 

static PyMethodDef PyTree_6dL_methods[] = {

  { (char*) "size",(PyCFunction)PyTree_6dL_size, METH_NOARGS, (char*) ""},

  { (char*) "optimize",(PyCFunction)PyTree_6dL_optimize, METH_NOARGS, (char*) ""},

  { (char*) "get_all",(PyCFunction)PyTree_6dL_get_all, METH_NOARGS, (char*) ""},

  { (char*) "add",(PyCFunction)PyTree_6dL_add, METH_VARARGS, (char*) ""},

  { (char*) "find_exact",(PyCFunction)PyTree_6dL_find_exact, METH_VARARGS, (char*) ""},

  { (char*) "find_nearest",(PyCFunction)PyTree_6dL_find_nearest, METH_VARARGS, (char*) ""},

  { (char*) "find_within_range",(PyCFunction)PyTree_6dL_find_within_range, METH_VARARGS, (char*) ""},

  { (char*) "count_within_range",(PyCFunction)PyTree_6dL_count_within_range, METH_VARARGS, (char*) ""},

  { (char*) "remove",(PyCFunction)PyTree_6dL_remove, METH_VARARGS, (char*) ""},

  {NULL, NULL, 0, NULL}  // Sentinel 
};


static PyTypeObject PyTree_6dLType = {
    PyObject_HEAD_INIT(NULL)
    0,                                // ob_size
    (char*)"KDTree_6Double",           // tp_name, __name__
    sizeof(PyTree_6dL),           // tp_basicsize
    0,                                // tp_itemsize
    (destructor)PyTree_6dL_dealloc, // tp_dealloc
    0,                                // tp_print
    0,                                // tp_getattr, __getattr__
    0,                                // tp_setattr, __setattr__
    0,                                // tp_compare
    (reprfunc)PyTree_6dL_repr,    // tp_repr, __repr__
    0,                                // tp_as_number
    0,                                // tp_as_sequence
    0,                                // tp_as_mapping
    0,                                // tp_hash
    0,                                // tp_call, __call__
    (reprfunc)PyTree_6dL_str,     // tp_str, __str__
    0,                                // tp_getattro
    0,                                // tp_setattro
    0,                                // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags
    (char*)"PyTree_6dL objects",  // tp_doc, __doc__
    0,		                      // tp_traverse 
    0,		                      // tp_clear 
    0,		                      // tp_richcompare 
    0,		                      // tp_weaklistoffset 
    0,		                      // tp_iter, __iter__
    0,		                      // tp_iternext 
    PyTree_6dL_methods,           // tp_methods 
    PyTree_6dL_members,           // tp_members 
    0,                                // tp_getset 
    0,                                // tp_base 
    0,                                // tp_dict, __dict__
    0,                                // tp_descr_get 
    0,                                // tp_descr_set 
    0,                                // tp_dictoffset 
    (initproc)PyTree_6dL_init,    // tp_init, __init__
    0,                                // tp_alloc 
    PyTree_6dL_new,               // tp_new 
};



#define RECORD_7dL record_t<7, double, unsigned long long >
#define RECORD_7dL_point record_t<7, double, unsigned long long >::point_t
#define RECORD_7dL_coord record_t<7, double, unsigned long long >::coord_t
#define KDTREE_TYPE_7dL KDTree::KDTree<7, RECORD_7dL, std::pointer_to_binary_function<RECORD_7dL,int,double> >
#define PYKDTREE_7dL PyKDTree<7, double, unsigned long long >




inline bool operator==(RECORD_7dL const& A, RECORD_7dL const& B) {
    return A.data == B.data && A.point[0] == B.point[0] && A.point[1] == B.point[1] && A.point[2] == B.point[2] && A.point[3] == B.point[3] && A.point[4] == B.point[4] && A.point[5] == B.point[5] && A.point[6] == B.point[6] ;
}

std::ostream& operator<<(std::ostream& out, RECORD_7dL const& T)
{
    return out << '(' << T.point[0] << ',' << T.point[1] << ',' << T.point[2] << ',' << T.point[3] << ',' << T.point[4] << ',' << T.point[5] << ',' << T.point[6] << '|' << T.data << ')';
}



// PyTree_7dL
// PYKDTREE_7dL
// RECORD_7dL
// RECORD_7dL_point
// RECORD_7dL_coord


typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  PYKDTREE_7dL *tree;
} PyTree_7dL;


static PyObject* PyTree_7dL_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyTree_7dL *self;
  
  self = (PyTree_7dL*)type->tp_alloc(type, 0);
  if(self)
    self->tree = NULL;
  
  return (PyObject*)self;
}


static int PyTree_7dL_init(PyTree_7dL *self, PyObject *args, PyObject *kwds)
{
  self->tree = (PYKDTREE_7dL *) new PYKDTREE_7dL();
  return 0;
}


// New datatype dealloc method
static void PyTree_7dL_dealloc(PyTree_7dL* self) {
  if(self->tree)
    delete self->tree;
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PyTree_7dL_repr(PyTree_7dL* self) {

  return PyString_FromString( (char*) "<KDTree_7Double instance>" );
}

static PyObject* PyTree_7dL_str(PyTree_7dL* self) {

  std::ostringstream output;
  std::vector<RECORD_7dL> *temp = 0;
  size_t size;
  size_t full = 6, half = 3;

  if(self && self->tree) {

    size = (int) self->tree->size();

    temp = (std::vector<RECORD_7dL > *) self->tree->get_all();
    
    if (size > full) {
      
      std::vector<RECORD_7dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < half; i++, iter++)
	output << (*iter) << std::endl;
      
      output << "..." << std::endl;
      
      for (size_t i = size-half-1; i < size; i++) 
	output << temp->at(i) << std::endl;
      
    }

    else if( size <= full ) {
    
      std::vector<RECORD_7dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < full; i++, iter++)
        output << (*iter) << std::endl;
    
    }
  
    else {

      output << "[]";
    }
    
    // delete the pointer after use!
    delete temp;
  }

  return PyString_FromString( output.str().c_str() );
}



static PyObject * PyTree_7dL_add(PyTree_7dL* self, PyObject *args, PyObject *kwds) {

  RECORD_7dL record;
  PYKDTREE_7dL *tree;

  if ( !PyArg_ParseTuple(args, "((ddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (7 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && ((PyTree_7dL*)self)->tree) {

    tree = (PYKDTREE_7dL *)((PyTree_7dL*)self)->tree;
    tree->add(record);

  } else {
    PyErr_SetString(PyExc_TypeError,"Adding record failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject * PyTree_7dL_remove(PyTree_7dL* self, PyObject *args, PyObject *kwds) {

  RECORD_7dL record;
  bool success;
  PyObject *result = Py_False;

  if ( !PyArg_ParseTuple(args, "((ddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (7 dim double vector, unsigned long long value)");
    return NULL;
  }

  if(self && self->tree) {
    success = self->tree->remove(record);
  } else {
    PyErr_SetString(PyExc_TypeError,"Removing record failed!");
    return NULL;
  }

  if(success)
    result = Py_True;

  Py_INCREF(result);
  return result;  
}

static PyObject * PyTree_7dL_size(PyTree_7dL* self) {

  int result;

  if(self && self->tree) {
    result = (int) self->tree->size();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing size() failed!");
    return NULL;
  }

  return Py_BuildValue((char*)"i", result);
}

static PyObject * PyTree_7dL_optimize(PyTree_7dL* self) {

  if(self && self->tree) {
    self->tree->optimize();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject * PyTree_7dL_find_exact(PyTree_7dL* self, PyObject *args, PyObject *kwds) {

  RECORD_7dL record;
  RECORD_7dL *temp = NULL;
  PyObject *result = NULL;

  if ( !PyArg_ParseTuple(args, "((ddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (7 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && self->tree) {

    temp = (RECORD_7dL *) self->tree->find_exact(record);

  } else {

    PyErr_SetString(PyExc_TypeError,"find exact failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(ddddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4],temp->point[5],temp->point[6]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}



static PyObject * PyTree_7dL_find_nearest(PyTree_7dL* self, PyObject *args, PyObject *kwds) {

  RECORD_7dL_point point_t;
  RECORD_7dL_coord *coord_t = NULL;
  RECORD_7dL *temp = NULL;
  PyObject *result = NULL;
  PYKDTREE_7dL *tree;

  if ( !PyArg_ParseTuple(args, "(ddddddd)",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6]
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 7 elements: 7 dim double vector");
    return NULL;
  }
  
  if(self && ((PyTree_7dL*)self)->tree) {

    tree = (PYKDTREE_7dL *)((PyTree_7dL*)self)->tree;

    coord_t = point_t;

    temp = (RECORD_7dL *) tree->find_nearest(coord_t);

  } else {

    PyErr_SetString(PyExc_TypeError,"find nearest failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(ddddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4],temp->point[5],temp->point[6]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}

static PyObject * PyTree_7dL_count_within_range(PyTree_7dL* self, PyObject *args, PyObject *kwds) {

  RECORD_7dL_point point_t;
  RECORD_7dL_coord *coord_t = NULL;
  size_t temp;
  PyObject *result = NULL;
  PYKDTREE_7dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(ddddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (7 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_7dL*)self)->tree) {

    tree = (PYKDTREE_7dL *)((PyTree_7dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (size_t) tree->count_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"count within range failed!");
    return NULL;
  }

  result = Py_BuildValue("i", (int)temp);
  return result;
}



static PyObject * PyTree_7dL_find_within_range(PyTree_7dL* self, PyObject *args, PyObject *kwds) {

  RECORD_7dL_point point_t;
  RECORD_7dL_coord *coord_t = NULL;
  std::vector<RECORD_7dL> *temp = 0;
  PYKDTREE_7dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(ddddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (7 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_7dL*)self)->tree) {

    tree = (PYKDTREE_7dL *)((PyTree_7dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (std::vector<RECORD_7dL > *) tree->find_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"find within range failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_7dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(ddddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4],(*iter).point[5],(*iter).point[6], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;
  
  return result;
}


static PyObject * PyTree_7dL_get_all(PyTree_7dL* self) {

  std::vector<RECORD_7dL> *temp = 0;

  if(self && self->tree) {

    temp = (std::vector<RECORD_7dL > *) self->tree->get_all();

  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_7dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(ddddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4],(*iter).point[5],(*iter).point[6], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;

  return result;
}




// Attributes of class are inserted and 'named' here  

static PyMemberDef PyTree_7dL_members[] = {
  // {(char*)"name", T_OBJECT_EX, offsetof(PyTree_7dL, name), 0,
  //  (char*)"name"},

  // {(char*)"datatype", T_INT, offsetof(PyTree_7dL, datatype), 0,
  //  (char*)"datatype"},

  // {(char*)"unit", T_OBJECT_EX, offsetof(PyTree_7dL, unit), 0,
  //  (char*)"unit"},

  // {(char*)"desc", T_OBJECT_EX, offsetof(PyTree_7dL, desc), 0,
  //  (char*)"desc"},
  {NULL}  // Sentinel 
};

// Methods of class are inserted and 'named' here  

// Notice that the function signature type is mapped with a flag 
// METH_NOARGS, METH_VARARGS, or METH_VARARGS | METH_KEYWORDS 

static PyMethodDef PyTree_7dL_methods[] = {

  { (char*) "size",(PyCFunction)PyTree_7dL_size, METH_NOARGS, (char*) ""},

  { (char*) "optimize",(PyCFunction)PyTree_7dL_optimize, METH_NOARGS, (char*) ""},

  { (char*) "get_all",(PyCFunction)PyTree_7dL_get_all, METH_NOARGS, (char*) ""},

  { (char*) "add",(PyCFunction)PyTree_7dL_add, METH_VARARGS, (char*) ""},

  { (char*) "find_exact",(PyCFunction)PyTree_7dL_find_exact, METH_VARARGS, (char*) ""},

  { (char*) "find_nearest",(PyCFunction)PyTree_7dL_find_nearest, METH_VARARGS, (char*) ""},

  { (char*) "find_within_range",(PyCFunction)PyTree_7dL_find_within_range, METH_VARARGS, (char*) ""},

  { (char*) "count_within_range",(PyCFunction)PyTree_7dL_count_within_range, METH_VARARGS, (char*) ""},

  { (char*) "remove",(PyCFunction)PyTree_7dL_remove, METH_VARARGS, (char*) ""},

  {NULL, NULL, 0, NULL}  // Sentinel 
};


static PyTypeObject PyTree_7dLType = {
    PyObject_HEAD_INIT(NULL)
    0,                                // ob_size
    (char*)"KDTree_7Double",           // tp_name, __name__
    sizeof(PyTree_7dL),           // tp_basicsize
    0,                                // tp_itemsize
    (destructor)PyTree_7dL_dealloc, // tp_dealloc
    0,                                // tp_print
    0,                                // tp_getattr, __getattr__
    0,                                // tp_setattr, __setattr__
    0,                                // tp_compare
    (reprfunc)PyTree_7dL_repr,    // tp_repr, __repr__
    0,                                // tp_as_number
    0,                                // tp_as_sequence
    0,                                // tp_as_mapping
    0,                                // tp_hash
    0,                                // tp_call, __call__
    (reprfunc)PyTree_7dL_str,     // tp_str, __str__
    0,                                // tp_getattro
    0,                                // tp_setattro
    0,                                // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags
    (char*)"PyTree_7dL objects",  // tp_doc, __doc__
    0,		                      // tp_traverse 
    0,		                      // tp_clear 
    0,		                      // tp_richcompare 
    0,		                      // tp_weaklistoffset 
    0,		                      // tp_iter, __iter__
    0,		                      // tp_iternext 
    PyTree_7dL_methods,           // tp_methods 
    PyTree_7dL_members,           // tp_members 
    0,                                // tp_getset 
    0,                                // tp_base 
    0,                                // tp_dict, __dict__
    0,                                // tp_descr_get 
    0,                                // tp_descr_set 
    0,                                // tp_dictoffset 
    (initproc)PyTree_7dL_init,    // tp_init, __init__
    0,                                // tp_alloc 
    PyTree_7dL_new,               // tp_new 
};



#define RECORD_8dL record_t<8, double, unsigned long long >
#define RECORD_8dL_point record_t<8, double, unsigned long long >::point_t
#define RECORD_8dL_coord record_t<8, double, unsigned long long >::coord_t
#define KDTREE_TYPE_8dL KDTree::KDTree<8, RECORD_8dL, std::pointer_to_binary_function<RECORD_8dL,int,double> >
#define PYKDTREE_8dL PyKDTree<8, double, unsigned long long >




inline bool operator==(RECORD_8dL const& A, RECORD_8dL const& B) {
    return A.data == B.data && A.point[0] == B.point[0] && A.point[1] == B.point[1] && A.point[2] == B.point[2] && A.point[3] == B.point[3] && A.point[4] == B.point[4] && A.point[5] == B.point[5] && A.point[6] == B.point[6] && A.point[7] == B.point[7] ;
}

std::ostream& operator<<(std::ostream& out, RECORD_8dL const& T)
{
    return out << '(' << T.point[0] << ',' << T.point[1] << ',' << T.point[2] << ',' << T.point[3] << ',' << T.point[4] << ',' << T.point[5] << ',' << T.point[6] << ',' << T.point[7] << '|' << T.data << ')';
}



// PyTree_8dL
// PYKDTREE_8dL
// RECORD_8dL
// RECORD_8dL_point
// RECORD_8dL_coord


typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  PYKDTREE_8dL *tree;
} PyTree_8dL;


static PyObject* PyTree_8dL_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyTree_8dL *self;
  
  self = (PyTree_8dL*)type->tp_alloc(type, 0);
  if(self)
    self->tree = NULL;
  
  return (PyObject*)self;
}


static int PyTree_8dL_init(PyTree_8dL *self, PyObject *args, PyObject *kwds)
{
  self->tree = (PYKDTREE_8dL *) new PYKDTREE_8dL();
  return 0;
}


// New datatype dealloc method
static void PyTree_8dL_dealloc(PyTree_8dL* self) {
  if(self->tree)
    delete self->tree;
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PyTree_8dL_repr(PyTree_8dL* self) {

  return PyString_FromString( (char*) "<KDTree_8Double instance>" );
}

static PyObject* PyTree_8dL_str(PyTree_8dL* self) {

  std::ostringstream output;
  std::vector<RECORD_8dL> *temp = 0;
  size_t size;
  size_t full = 6, half = 3;

  if(self && self->tree) {

    size = (int) self->tree->size();

    temp = (std::vector<RECORD_8dL > *) self->tree->get_all();
    
    if (size > full) {
      
      std::vector<RECORD_8dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < half; i++, iter++)
	output << (*iter) << std::endl;
      
      output << "..." << std::endl;
      
      for (size_t i = size-half-1; i < size; i++) 
	output << temp->at(i) << std::endl;
      
    }

    else if( size <= full ) {
    
      std::vector<RECORD_8dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < full; i++, iter++)
        output << (*iter) << std::endl;
    
    }
  
    else {

      output << "[]";
    }
    
    // delete the pointer after use!
    delete temp;
  }

  return PyString_FromString( output.str().c_str() );
}



static PyObject * PyTree_8dL_add(PyTree_8dL* self, PyObject *args, PyObject *kwds) {

  RECORD_8dL record;
  PYKDTREE_8dL *tree;

  if ( !PyArg_ParseTuple(args, "((dddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
&record.point[7],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (8 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && ((PyTree_8dL*)self)->tree) {

    tree = (PYKDTREE_8dL *)((PyTree_8dL*)self)->tree;
    tree->add(record);

  } else {
    PyErr_SetString(PyExc_TypeError,"Adding record failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject * PyTree_8dL_remove(PyTree_8dL* self, PyObject *args, PyObject *kwds) {

  RECORD_8dL record;
  bool success;
  PyObject *result = Py_False;

  if ( !PyArg_ParseTuple(args, "((dddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
&record.point[7],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (8 dim double vector, unsigned long long value)");
    return NULL;
  }

  if(self && self->tree) {
    success = self->tree->remove(record);
  } else {
    PyErr_SetString(PyExc_TypeError,"Removing record failed!");
    return NULL;
  }

  if(success)
    result = Py_True;

  Py_INCREF(result);
  return result;  
}

static PyObject * PyTree_8dL_size(PyTree_8dL* self) {

  int result;

  if(self && self->tree) {
    result = (int) self->tree->size();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing size() failed!");
    return NULL;
  }

  return Py_BuildValue((char*)"i", result);
}

static PyObject * PyTree_8dL_optimize(PyTree_8dL* self) {

  if(self && self->tree) {
    self->tree->optimize();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject * PyTree_8dL_find_exact(PyTree_8dL* self, PyObject *args, PyObject *kwds) {

  RECORD_8dL record;
  RECORD_8dL *temp = NULL;
  PyObject *result = NULL;

  if ( !PyArg_ParseTuple(args, "((dddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
&record.point[7],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (8 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && self->tree) {

    temp = (RECORD_8dL *) self->tree->find_exact(record);

  } else {

    PyErr_SetString(PyExc_TypeError,"find exact failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(dddddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4],temp->point[5],temp->point[6],temp->point[7]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}



static PyObject * PyTree_8dL_find_nearest(PyTree_8dL* self, PyObject *args, PyObject *kwds) {

  RECORD_8dL_point point_t;
  RECORD_8dL_coord *coord_t = NULL;
  RECORD_8dL *temp = NULL;
  PyObject *result = NULL;
  PYKDTREE_8dL *tree;

  if ( !PyArg_ParseTuple(args, "(dddddddd)",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
&point_t[7]
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 8 elements: 8 dim double vector");
    return NULL;
  }
  
  if(self && ((PyTree_8dL*)self)->tree) {

    tree = (PYKDTREE_8dL *)((PyTree_8dL*)self)->tree;

    coord_t = point_t;

    temp = (RECORD_8dL *) tree->find_nearest(coord_t);

  } else {

    PyErr_SetString(PyExc_TypeError,"find nearest failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(dddddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4],temp->point[5],temp->point[6],temp->point[7]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}

static PyObject * PyTree_8dL_count_within_range(PyTree_8dL* self, PyObject *args, PyObject *kwds) {

  RECORD_8dL_point point_t;
  RECORD_8dL_coord *coord_t = NULL;
  size_t temp;
  PyObject *result = NULL;
  PYKDTREE_8dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(dddddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
&point_t[7],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (8 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_8dL*)self)->tree) {

    tree = (PYKDTREE_8dL *)((PyTree_8dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (size_t) tree->count_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"count within range failed!");
    return NULL;
  }

  result = Py_BuildValue("i", (int)temp);
  return result;
}



static PyObject * PyTree_8dL_find_within_range(PyTree_8dL* self, PyObject *args, PyObject *kwds) {

  RECORD_8dL_point point_t;
  RECORD_8dL_coord *coord_t = NULL;
  std::vector<RECORD_8dL> *temp = 0;
  PYKDTREE_8dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(dddddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
&point_t[7],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (8 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_8dL*)self)->tree) {

    tree = (PYKDTREE_8dL *)((PyTree_8dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (std::vector<RECORD_8dL > *) tree->find_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"find within range failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_8dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(dddddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4],(*iter).point[5],(*iter).point[6],(*iter).point[7], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;
  
  return result;
}


static PyObject * PyTree_8dL_get_all(PyTree_8dL* self) {

  std::vector<RECORD_8dL> *temp = 0;

  if(self && self->tree) {

    temp = (std::vector<RECORD_8dL > *) self->tree->get_all();

  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_8dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(dddddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4],(*iter).point[5],(*iter).point[6],(*iter).point[7], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;

  return result;
}




// Attributes of class are inserted and 'named' here  

static PyMemberDef PyTree_8dL_members[] = {
  // {(char*)"name", T_OBJECT_EX, offsetof(PyTree_8dL, name), 0,
  //  (char*)"name"},

  // {(char*)"datatype", T_INT, offsetof(PyTree_8dL, datatype), 0,
  //  (char*)"datatype"},

  // {(char*)"unit", T_OBJECT_EX, offsetof(PyTree_8dL, unit), 0,
  //  (char*)"unit"},

  // {(char*)"desc", T_OBJECT_EX, offsetof(PyTree_8dL, desc), 0,
  //  (char*)"desc"},
  {NULL}  // Sentinel 
};

// Methods of class are inserted and 'named' here  

// Notice that the function signature type is mapped with a flag 
// METH_NOARGS, METH_VARARGS, or METH_VARARGS | METH_KEYWORDS 

static PyMethodDef PyTree_8dL_methods[] = {

  { (char*) "size",(PyCFunction)PyTree_8dL_size, METH_NOARGS, (char*) ""},

  { (char*) "optimize",(PyCFunction)PyTree_8dL_optimize, METH_NOARGS, (char*) ""},

  { (char*) "get_all",(PyCFunction)PyTree_8dL_get_all, METH_NOARGS, (char*) ""},

  { (char*) "add",(PyCFunction)PyTree_8dL_add, METH_VARARGS, (char*) ""},

  { (char*) "find_exact",(PyCFunction)PyTree_8dL_find_exact, METH_VARARGS, (char*) ""},

  { (char*) "find_nearest",(PyCFunction)PyTree_8dL_find_nearest, METH_VARARGS, (char*) ""},

  { (char*) "find_within_range",(PyCFunction)PyTree_8dL_find_within_range, METH_VARARGS, (char*) ""},

  { (char*) "count_within_range",(PyCFunction)PyTree_8dL_count_within_range, METH_VARARGS, (char*) ""},

  { (char*) "remove",(PyCFunction)PyTree_8dL_remove, METH_VARARGS, (char*) ""},

  {NULL, NULL, 0, NULL}  // Sentinel 
};


static PyTypeObject PyTree_8dLType = {
    PyObject_HEAD_INIT(NULL)
    0,                                // ob_size
    (char*)"KDTree_8Double",           // tp_name, __name__
    sizeof(PyTree_8dL),           // tp_basicsize
    0,                                // tp_itemsize
    (destructor)PyTree_8dL_dealloc, // tp_dealloc
    0,                                // tp_print
    0,                                // tp_getattr, __getattr__
    0,                                // tp_setattr, __setattr__
    0,                                // tp_compare
    (reprfunc)PyTree_8dL_repr,    // tp_repr, __repr__
    0,                                // tp_as_number
    0,                                // tp_as_sequence
    0,                                // tp_as_mapping
    0,                                // tp_hash
    0,                                // tp_call, __call__
    (reprfunc)PyTree_8dL_str,     // tp_str, __str__
    0,                                // tp_getattro
    0,                                // tp_setattro
    0,                                // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags
    (char*)"PyTree_8dL objects",  // tp_doc, __doc__
    0,		                      // tp_traverse 
    0,		                      // tp_clear 
    0,		                      // tp_richcompare 
    0,		                      // tp_weaklistoffset 
    0,		                      // tp_iter, __iter__
    0,		                      // tp_iternext 
    PyTree_8dL_methods,           // tp_methods 
    PyTree_8dL_members,           // tp_members 
    0,                                // tp_getset 
    0,                                // tp_base 
    0,                                // tp_dict, __dict__
    0,                                // tp_descr_get 
    0,                                // tp_descr_set 
    0,                                // tp_dictoffset 
    (initproc)PyTree_8dL_init,    // tp_init, __init__
    0,                                // tp_alloc 
    PyTree_8dL_new,               // tp_new 
};



#define RECORD_9dL record_t<9, double, unsigned long long >
#define RECORD_9dL_point record_t<9, double, unsigned long long >::point_t
#define RECORD_9dL_coord record_t<9, double, unsigned long long >::coord_t
#define KDTREE_TYPE_9dL KDTree::KDTree<9, RECORD_9dL, std::pointer_to_binary_function<RECORD_9dL,int,double> >
#define PYKDTREE_9dL PyKDTree<9, double, unsigned long long >




inline bool operator==(RECORD_9dL const& A, RECORD_9dL const& B) {
    return A.data == B.data && A.point[0] == B.point[0] && A.point[1] == B.point[1] && A.point[2] == B.point[2] && A.point[3] == B.point[3] && A.point[4] == B.point[4] && A.point[5] == B.point[5] && A.point[6] == B.point[6] && A.point[7] == B.point[7] && A.point[8] == B.point[8] ;
}

std::ostream& operator<<(std::ostream& out, RECORD_9dL const& T)
{
    return out << '(' << T.point[0] << ',' << T.point[1] << ',' << T.point[2] << ',' << T.point[3] << ',' << T.point[4] << ',' << T.point[5] << ',' << T.point[6] << ',' << T.point[7] << ',' << T.point[8] << '|' << T.data << ')';
}



// PyTree_9dL
// PYKDTREE_9dL
// RECORD_9dL
// RECORD_9dL_point
// RECORD_9dL_coord


typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  PYKDTREE_9dL *tree;
} PyTree_9dL;


static PyObject* PyTree_9dL_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyTree_9dL *self;
  
  self = (PyTree_9dL*)type->tp_alloc(type, 0);
  if(self)
    self->tree = NULL;
  
  return (PyObject*)self;
}


static int PyTree_9dL_init(PyTree_9dL *self, PyObject *args, PyObject *kwds)
{
  self->tree = (PYKDTREE_9dL *) new PYKDTREE_9dL();
  return 0;
}


// New datatype dealloc method
static void PyTree_9dL_dealloc(PyTree_9dL* self) {
  if(self->tree)
    delete self->tree;
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PyTree_9dL_repr(PyTree_9dL* self) {

  return PyString_FromString( (char*) "<KDTree_9Double instance>" );
}

static PyObject* PyTree_9dL_str(PyTree_9dL* self) {

  std::ostringstream output;
  std::vector<RECORD_9dL> *temp = 0;
  size_t size;
  size_t full = 6, half = 3;

  if(self && self->tree) {

    size = (int) self->tree->size();

    temp = (std::vector<RECORD_9dL > *) self->tree->get_all();
    
    if (size > full) {
      
      std::vector<RECORD_9dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < half; i++, iter++)
	output << (*iter) << std::endl;
      
      output << "..." << std::endl;
      
      for (size_t i = size-half-1; i < size; i++) 
	output << temp->at(i) << std::endl;
      
    }

    else if( size <= full ) {
    
      std::vector<RECORD_9dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < full; i++, iter++)
        output << (*iter) << std::endl;
    
    }
  
    else {

      output << "[]";
    }
    
    // delete the pointer after use!
    delete temp;
  }

  return PyString_FromString( output.str().c_str() );
}



static PyObject * PyTree_9dL_add(PyTree_9dL* self, PyObject *args, PyObject *kwds) {

  RECORD_9dL record;
  PYKDTREE_9dL *tree;

  if ( !PyArg_ParseTuple(args, "((ddddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
&record.point[7],
&record.point[8],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (9 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && ((PyTree_9dL*)self)->tree) {

    tree = (PYKDTREE_9dL *)((PyTree_9dL*)self)->tree;
    tree->add(record);

  } else {
    PyErr_SetString(PyExc_TypeError,"Adding record failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject * PyTree_9dL_remove(PyTree_9dL* self, PyObject *args, PyObject *kwds) {

  RECORD_9dL record;
  bool success;
  PyObject *result = Py_False;

  if ( !PyArg_ParseTuple(args, "((ddddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
&record.point[7],
&record.point[8],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (9 dim double vector, unsigned long long value)");
    return NULL;
  }

  if(self && self->tree) {
    success = self->tree->remove(record);
  } else {
    PyErr_SetString(PyExc_TypeError,"Removing record failed!");
    return NULL;
  }

  if(success)
    result = Py_True;

  Py_INCREF(result);
  return result;  
}

static PyObject * PyTree_9dL_size(PyTree_9dL* self) {

  int result;

  if(self && self->tree) {
    result = (int) self->tree->size();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing size() failed!");
    return NULL;
  }

  return Py_BuildValue((char*)"i", result);
}

static PyObject * PyTree_9dL_optimize(PyTree_9dL* self) {

  if(self && self->tree) {
    self->tree->optimize();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject * PyTree_9dL_find_exact(PyTree_9dL* self, PyObject *args, PyObject *kwds) {

  RECORD_9dL record;
  RECORD_9dL *temp = NULL;
  PyObject *result = NULL;

  if ( !PyArg_ParseTuple(args, "((ddddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
&record.point[7],
&record.point[8],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (9 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && self->tree) {

    temp = (RECORD_9dL *) self->tree->find_exact(record);

  } else {

    PyErr_SetString(PyExc_TypeError,"find exact failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(ddddddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4],temp->point[5],temp->point[6],temp->point[7],temp->point[8]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}



static PyObject * PyTree_9dL_find_nearest(PyTree_9dL* self, PyObject *args, PyObject *kwds) {

  RECORD_9dL_point point_t;
  RECORD_9dL_coord *coord_t = NULL;
  RECORD_9dL *temp = NULL;
  PyObject *result = NULL;
  PYKDTREE_9dL *tree;

  if ( !PyArg_ParseTuple(args, "(ddddddddd)",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
&point_t[7],
&point_t[8]
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 9 elements: 9 dim double vector");
    return NULL;
  }
  
  if(self && ((PyTree_9dL*)self)->tree) {

    tree = (PYKDTREE_9dL *)((PyTree_9dL*)self)->tree;

    coord_t = point_t;

    temp = (RECORD_9dL *) tree->find_nearest(coord_t);

  } else {

    PyErr_SetString(PyExc_TypeError,"find nearest failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(ddddddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4],temp->point[5],temp->point[6],temp->point[7],temp->point[8]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}

static PyObject * PyTree_9dL_count_within_range(PyTree_9dL* self, PyObject *args, PyObject *kwds) {

  RECORD_9dL_point point_t;
  RECORD_9dL_coord *coord_t = NULL;
  size_t temp;
  PyObject *result = NULL;
  PYKDTREE_9dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(ddddddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
&point_t[7],
&point_t[8],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (9 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_9dL*)self)->tree) {

    tree = (PYKDTREE_9dL *)((PyTree_9dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (size_t) tree->count_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"count within range failed!");
    return NULL;
  }

  result = Py_BuildValue("i", (int)temp);
  return result;
}



static PyObject * PyTree_9dL_find_within_range(PyTree_9dL* self, PyObject *args, PyObject *kwds) {

  RECORD_9dL_point point_t;
  RECORD_9dL_coord *coord_t = NULL;
  std::vector<RECORD_9dL> *temp = 0;
  PYKDTREE_9dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(ddddddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
&point_t[7],
&point_t[8],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (9 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_9dL*)self)->tree) {

    tree = (PYKDTREE_9dL *)((PyTree_9dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (std::vector<RECORD_9dL > *) tree->find_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"find within range failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_9dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(ddddddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4],(*iter).point[5],(*iter).point[6],(*iter).point[7],(*iter).point[8], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;
  
  return result;
}


static PyObject * PyTree_9dL_get_all(PyTree_9dL* self) {

  std::vector<RECORD_9dL> *temp = 0;

  if(self && self->tree) {

    temp = (std::vector<RECORD_9dL > *) self->tree->get_all();

  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_9dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(ddddddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4],(*iter).point[5],(*iter).point[6],(*iter).point[7],(*iter).point[8], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;

  return result;
}




// Attributes of class are inserted and 'named' here  

static PyMemberDef PyTree_9dL_members[] = {
  // {(char*)"name", T_OBJECT_EX, offsetof(PyTree_9dL, name), 0,
  //  (char*)"name"},

  // {(char*)"datatype", T_INT, offsetof(PyTree_9dL, datatype), 0,
  //  (char*)"datatype"},

  // {(char*)"unit", T_OBJECT_EX, offsetof(PyTree_9dL, unit), 0,
  //  (char*)"unit"},

  // {(char*)"desc", T_OBJECT_EX, offsetof(PyTree_9dL, desc), 0,
  //  (char*)"desc"},
  {NULL}  // Sentinel 
};

// Methods of class are inserted and 'named' here  

// Notice that the function signature type is mapped with a flag 
// METH_NOARGS, METH_VARARGS, or METH_VARARGS | METH_KEYWORDS 

static PyMethodDef PyTree_9dL_methods[] = {

  { (char*) "size",(PyCFunction)PyTree_9dL_size, METH_NOARGS, (char*) ""},

  { (char*) "optimize",(PyCFunction)PyTree_9dL_optimize, METH_NOARGS, (char*) ""},

  { (char*) "get_all",(PyCFunction)PyTree_9dL_get_all, METH_NOARGS, (char*) ""},

  { (char*) "add",(PyCFunction)PyTree_9dL_add, METH_VARARGS, (char*) ""},

  { (char*) "find_exact",(PyCFunction)PyTree_9dL_find_exact, METH_VARARGS, (char*) ""},

  { (char*) "find_nearest",(PyCFunction)PyTree_9dL_find_nearest, METH_VARARGS, (char*) ""},

  { (char*) "find_within_range",(PyCFunction)PyTree_9dL_find_within_range, METH_VARARGS, (char*) ""},

  { (char*) "count_within_range",(PyCFunction)PyTree_9dL_count_within_range, METH_VARARGS, (char*) ""},

  { (char*) "remove",(PyCFunction)PyTree_9dL_remove, METH_VARARGS, (char*) ""},

  {NULL, NULL, 0, NULL}  // Sentinel 
};


static PyTypeObject PyTree_9dLType = {
    PyObject_HEAD_INIT(NULL)
    0,                                // ob_size
    (char*)"KDTree_9Double",           // tp_name, __name__
    sizeof(PyTree_9dL),           // tp_basicsize
    0,                                // tp_itemsize
    (destructor)PyTree_9dL_dealloc, // tp_dealloc
    0,                                // tp_print
    0,                                // tp_getattr, __getattr__
    0,                                // tp_setattr, __setattr__
    0,                                // tp_compare
    (reprfunc)PyTree_9dL_repr,    // tp_repr, __repr__
    0,                                // tp_as_number
    0,                                // tp_as_sequence
    0,                                // tp_as_mapping
    0,                                // tp_hash
    0,                                // tp_call, __call__
    (reprfunc)PyTree_9dL_str,     // tp_str, __str__
    0,                                // tp_getattro
    0,                                // tp_setattro
    0,                                // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags
    (char*)"PyTree_9dL objects",  // tp_doc, __doc__
    0,		                      // tp_traverse 
    0,		                      // tp_clear 
    0,		                      // tp_richcompare 
    0,		                      // tp_weaklistoffset 
    0,		                      // tp_iter, __iter__
    0,		                      // tp_iternext 
    PyTree_9dL_methods,           // tp_methods 
    PyTree_9dL_members,           // tp_members 
    0,                                // tp_getset 
    0,                                // tp_base 
    0,                                // tp_dict, __dict__
    0,                                // tp_descr_get 
    0,                                // tp_descr_set 
    0,                                // tp_dictoffset 
    (initproc)PyTree_9dL_init,    // tp_init, __init__
    0,                                // tp_alloc 
    PyTree_9dL_new,               // tp_new 
};



#define RECORD_10dL record_t<10, double, unsigned long long >
#define RECORD_10dL_point record_t<10, double, unsigned long long >::point_t
#define RECORD_10dL_coord record_t<10, double, unsigned long long >::coord_t
#define KDTREE_TYPE_10dL KDTree::KDTree<10, RECORD_10dL, std::pointer_to_binary_function<RECORD_10dL,int,double> >
#define PYKDTREE_10dL PyKDTree<10, double, unsigned long long >




inline bool operator==(RECORD_10dL const& A, RECORD_10dL const& B) {
    return A.data == B.data && A.point[0] == B.point[0] && A.point[1] == B.point[1] && A.point[2] == B.point[2] && A.point[3] == B.point[3] && A.point[4] == B.point[4] && A.point[5] == B.point[5] && A.point[6] == B.point[6] && A.point[7] == B.point[7] && A.point[8] == B.point[8] && A.point[9] == B.point[9] ;
}

std::ostream& operator<<(std::ostream& out, RECORD_10dL const& T)
{
    return out << '(' << T.point[0] << ',' << T.point[1] << ',' << T.point[2] << ',' << T.point[3] << ',' << T.point[4] << ',' << T.point[5] << ',' << T.point[6] << ',' << T.point[7] << ',' << T.point[8] << ',' << T.point[9] << '|' << T.data << ')';
}



// PyTree_10dL
// PYKDTREE_10dL
// RECORD_10dL
// RECORD_10dL_point
// RECORD_10dL_coord


typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  PYKDTREE_10dL *tree;
} PyTree_10dL;


static PyObject* PyTree_10dL_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyTree_10dL *self;
  
  self = (PyTree_10dL*)type->tp_alloc(type, 0);
  if(self)
    self->tree = NULL;
  
  return (PyObject*)self;
}


static int PyTree_10dL_init(PyTree_10dL *self, PyObject *args, PyObject *kwds)
{
  self->tree = (PYKDTREE_10dL *) new PYKDTREE_10dL();
  return 0;
}


// New datatype dealloc method
static void PyTree_10dL_dealloc(PyTree_10dL* self) {
  if(self->tree)
    delete self->tree;
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* PyTree_10dL_repr(PyTree_10dL* self) {

  return PyString_FromString( (char*) "<KDTree_10Double instance>" );
}

static PyObject* PyTree_10dL_str(PyTree_10dL* self) {

  std::ostringstream output;
  std::vector<RECORD_10dL> *temp = 0;
  size_t size;
  size_t full = 6, half = 3;

  if(self && self->tree) {

    size = (int) self->tree->size();

    temp = (std::vector<RECORD_10dL > *) self->tree->get_all();
    
    if (size > full) {
      
      std::vector<RECORD_10dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < half; i++, iter++)
	output << (*iter) << std::endl;
      
      output << "..." << std::endl;
      
      for (size_t i = size-half-1; i < size; i++) 
	output << temp->at(i) << std::endl;
      
    }

    else if( size <= full ) {
    
      std::vector<RECORD_10dL >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < full; i++, iter++)
        output << (*iter) << std::endl;
    
    }
  
    else {

      output << "[]";
    }
    
    // delete the pointer after use!
    delete temp;
  }

  return PyString_FromString( output.str().c_str() );
}



static PyObject * PyTree_10dL_add(PyTree_10dL* self, PyObject *args, PyObject *kwds) {

  RECORD_10dL record;
  PYKDTREE_10dL *tree;

  if ( !PyArg_ParseTuple(args, "((dddddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
&record.point[7],
&record.point[8],
&record.point[9],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (10 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && ((PyTree_10dL*)self)->tree) {

    tree = (PYKDTREE_10dL *)((PyTree_10dL*)self)->tree;
    tree->add(record);

  } else {
    PyErr_SetString(PyExc_TypeError,"Adding record failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject * PyTree_10dL_remove(PyTree_10dL* self, PyObject *args, PyObject *kwds) {

  RECORD_10dL record;
  bool success;
  PyObject *result = Py_False;

  if ( !PyArg_ParseTuple(args, "((dddddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
&record.point[7],
&record.point[8],
&record.point[9],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (10 dim double vector, unsigned long long value)");
    return NULL;
  }

  if(self && self->tree) {
    success = self->tree->remove(record);
  } else {
    PyErr_SetString(PyExc_TypeError,"Removing record failed!");
    return NULL;
  }

  if(success)
    result = Py_True;

  Py_INCREF(result);
  return result;  
}

static PyObject * PyTree_10dL_size(PyTree_10dL* self) {

  int result;

  if(self && self->tree) {
    result = (int) self->tree->size();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing size() failed!");
    return NULL;
  }

  return Py_BuildValue((char*)"i", result);
}

static PyObject * PyTree_10dL_optimize(PyTree_10dL* self) {

  if(self && self->tree) {
    self->tree->optimize();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject * PyTree_10dL_find_exact(PyTree_10dL* self, PyObject *args, PyObject *kwds) {

  RECORD_10dL record;
  RECORD_10dL *temp = NULL;
  PyObject *result = NULL;

  if ( !PyArg_ParseTuple(args, "((dddddddddd)L)",
			 &record.point[0],
&record.point[1],
&record.point[2],
&record.point[3],
&record.point[4],
&record.point[5],
&record.point[6],
&record.point[7],
&record.point[8],
&record.point[9],
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (10 dim double vector, unsigned long long value)");
    return NULL;
  }
  
  if(self && self->tree) {

    temp = (RECORD_10dL *) self->tree->find_exact(record);

  } else {

    PyErr_SetString(PyExc_TypeError,"find exact failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(dddddddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4],temp->point[5],temp->point[6],temp->point[7],temp->point[8],temp->point[9]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}



static PyObject * PyTree_10dL_find_nearest(PyTree_10dL* self, PyObject *args, PyObject *kwds) {

  RECORD_10dL_point point_t;
  RECORD_10dL_coord *coord_t = NULL;
  RECORD_10dL *temp = NULL;
  PyObject *result = NULL;
  PYKDTREE_10dL *tree;

  if ( !PyArg_ParseTuple(args, "(dddddddddd)",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
&point_t[7],
&point_t[8],
&point_t[9]
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 10 elements: 10 dim double vector");
    return NULL;
  }
  
  if(self && ((PyTree_10dL*)self)->tree) {

    tree = (PYKDTREE_10dL *)((PyTree_10dL*)self)->tree;

    coord_t = point_t;

    temp = (RECORD_10dL *) tree->find_nearest(coord_t);

  } else {

    PyErr_SetString(PyExc_TypeError,"find nearest failed!");
    return NULL;
  }

  if (temp != NULL) {

    result = PyTuple_New(2);

    if (!result) {
      PyErr_SetString(PyErr_Occurred(),"unable to create a tuple.");
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(dddddddddd)",
						 temp->point[0],temp->point[1],temp->point[2],temp->point[3],temp->point[4],temp->point[5],temp->point[6],temp->point[7],temp->point[8],temp->point[9]))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("L", temp->data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(b) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }

    // delete the pointer after use!
    delete temp;

  } else {
    result = Py_BuildValue("");
  }
  
  return result;
}

static PyObject * PyTree_10dL_count_within_range(PyTree_10dL* self, PyObject *args, PyObject *kwds) {

  RECORD_10dL_point point_t;
  RECORD_10dL_coord *coord_t = NULL;
  size_t temp;
  PyObject *result = NULL;
  PYKDTREE_10dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(dddddddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
&point_t[7],
&point_t[8],
&point_t[9],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (10 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_10dL*)self)->tree) {

    tree = (PYKDTREE_10dL *)((PyTree_10dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (size_t) tree->count_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"count within range failed!");
    return NULL;
  }

  result = Py_BuildValue("i", (int)temp);
  return result;
}



static PyObject * PyTree_10dL_find_within_range(PyTree_10dL* self, PyObject *args, PyObject *kwds) {

  RECORD_10dL_point point_t;
  RECORD_10dL_coord *coord_t = NULL;
  std::vector<RECORD_10dL> *temp = 0;
  PYKDTREE_10dL *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(dddddddddd)d",
			 &point_t[0],
&point_t[1],
&point_t[2],
&point_t[3],
&point_t[4],
&point_t[5],
&point_t[6],
&point_t[7],
&point_t[8],
&point_t[9],
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (10 dim double vector, double value)");
    return NULL;
  }
  
  if(self && ((PyTree_10dL*)self)->tree) {

    tree = (PYKDTREE_10dL *)((PyTree_10dL*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (std::vector<RECORD_10dL > *) tree->find_within_range(coord_t, range);

  } else {

    PyErr_SetString(PyExc_TypeError,"find within range failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_10dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(dddddddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4],(*iter).point[5],(*iter).point[6],(*iter).point[7],(*iter).point[8],(*iter).point[9], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;
  
  return result;
}


static PyObject * PyTree_10dL_get_all(PyTree_10dL* self) {

  std::vector<RECORD_10dL> *temp = 0;

  if(self && self->tree) {

    temp = (std::vector<RECORD_10dL > *) self->tree->get_all();

  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  if(!temp) {
    PyErr_SetString(PyExc_TypeError,"function get_all() failed!");
    return NULL;
  }

  PyObject* result = PyList_New(temp->size());

  if (!result) {
    
    PyErr_SetString(PyErr_Occurred(),"unable to create a list.");
    return NULL;
  }

  std::vector<RECORD_10dL >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(dddddddddd)L",
	(*iter).point[0],(*iter).point[1],(*iter).point[2],(*iter).point[3],(*iter).point[4],(*iter).point[5],(*iter).point[6],(*iter).point[7],(*iter).point[8],(*iter).point[9], (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;

  return result;
}




// Attributes of class are inserted and 'named' here  

static PyMemberDef PyTree_10dL_members[] = {
  // {(char*)"name", T_OBJECT_EX, offsetof(PyTree_10dL, name), 0,
  //  (char*)"name"},

  // {(char*)"datatype", T_INT, offsetof(PyTree_10dL, datatype), 0,
  //  (char*)"datatype"},

  // {(char*)"unit", T_OBJECT_EX, offsetof(PyTree_10dL, unit), 0,
  //  (char*)"unit"},

  // {(char*)"desc", T_OBJECT_EX, offsetof(PyTree_10dL, desc), 0,
  //  (char*)"desc"},
  {NULL}  // Sentinel 
};

// Methods of class are inserted and 'named' here  

// Notice that the function signature type is mapped with a flag 
// METH_NOARGS, METH_VARARGS, or METH_VARARGS | METH_KEYWORDS 

static PyMethodDef PyTree_10dL_methods[] = {

  { (char*) "size",(PyCFunction)PyTree_10dL_size, METH_NOARGS, (char*) ""},

  { (char*) "optimize",(PyCFunction)PyTree_10dL_optimize, METH_NOARGS, (char*) ""},

  { (char*) "get_all",(PyCFunction)PyTree_10dL_get_all, METH_NOARGS, (char*) ""},

  { (char*) "add",(PyCFunction)PyTree_10dL_add, METH_VARARGS, (char*) ""},

  { (char*) "find_exact",(PyCFunction)PyTree_10dL_find_exact, METH_VARARGS, (char*) ""},

  { (char*) "find_nearest",(PyCFunction)PyTree_10dL_find_nearest, METH_VARARGS, (char*) ""},

  { (char*) "find_within_range",(PyCFunction)PyTree_10dL_find_within_range, METH_VARARGS, (char*) ""},

  { (char*) "count_within_range",(PyCFunction)PyTree_10dL_count_within_range, METH_VARARGS, (char*) ""},

  { (char*) "remove",(PyCFunction)PyTree_10dL_remove, METH_VARARGS, (char*) ""},

  {NULL, NULL, 0, NULL}  // Sentinel 
};


static PyTypeObject PyTree_10dLType = {
    PyObject_HEAD_INIT(NULL)
    0,                                // ob_size
    (char*)"KDTree_10Double",           // tp_name, __name__
    sizeof(PyTree_10dL),           // tp_basicsize
    0,                                // tp_itemsize
    (destructor)PyTree_10dL_dealloc, // tp_dealloc
    0,                                // tp_print
    0,                                // tp_getattr, __getattr__
    0,                                // tp_setattr, __setattr__
    0,                                // tp_compare
    (reprfunc)PyTree_10dL_repr,    // tp_repr, __repr__
    0,                                // tp_as_number
    0,                                // tp_as_sequence
    0,                                // tp_as_mapping
    0,                                // tp_hash
    0,                                // tp_call, __call__
    (reprfunc)PyTree_10dL_str,     // tp_str, __str__
    0,                                // tp_getattro
    0,                                // tp_setattro
    0,                                // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags
    (char*)"PyTree_10dL objects",  // tp_doc, __doc__
    0,		                      // tp_traverse 
    0,		                      // tp_clear 
    0,		                      // tp_richcompare 
    0,		                      // tp_weaklistoffset 
    0,		                      // tp_iter, __iter__
    0,		                      // tp_iternext 
    PyTree_10dL_methods,           // tp_methods 
    PyTree_10dL_members,           // tp_members 
    0,                                // tp_getset 
    0,                                // tp_base 
    0,                                // tp_dict, __dict__
    0,                                // tp_descr_get 
    0,                                // tp_descr_set 
    0,                                // tp_dictoffset 
    (initproc)PyTree_10dL_init,    // tp_init, __init__
    0,                                // tp_alloc 
    PyTree_10dL_new,               // tp_new 
};



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

  if (PyType_Ready(&PyTree_2dLType) < 0)
    return;

else if (PyType_Ready(&PyTree_3dLType) < 0)
    return;

else if (PyType_Ready(&PyTree_4dLType) < 0)
    return;

else if (PyType_Ready(&PyTree_5dLType) < 0)
    return;

else if (PyType_Ready(&PyTree_6dLType) < 0)
    return;

else if (PyType_Ready(&PyTree_7dLType) < 0)
    return;

else if (PyType_Ready(&PyTree_8dLType) < 0)
    return;

else if (PyType_Ready(&PyTree_9dLType) < 0)
    return;

else if (PyType_Ready(&PyTree_10dLType) < 0)
    return;


  //import_array();

  m = Py_InitModule3((char*)"_pykdtree", module_methods,
		     (char*)"Example module that creates an extension type.");
  
  if (m == NULL)
    return;
  
  Py_INCREF(&PyTree_2dLType);
Py_INCREF(&PyTree_3dLType);
Py_INCREF(&PyTree_4dLType);
Py_INCREF(&PyTree_5dLType);
Py_INCREF(&PyTree_6dLType);
Py_INCREF(&PyTree_7dLType);
Py_INCREF(&PyTree_8dLType);
Py_INCREF(&PyTree_9dLType);
Py_INCREF(&PyTree_10dLType);

  // C++ classes mapped to Python types go here and are then inserted into
  // the module namespace as Python classes
  PyModule_AddObject(m, (char*)"KDTree_2Double", (PyObject *)&PyTree_2dLType);
PyModule_AddObject(m, (char*)"KDTree_3Double", (PyObject *)&PyTree_3dLType);
PyModule_AddObject(m, (char*)"KDTree_4Double", (PyObject *)&PyTree_4dLType);
PyModule_AddObject(m, (char*)"KDTree_5Double", (PyObject *)&PyTree_5dLType);
PyModule_AddObject(m, (char*)"KDTree_6Double", (PyObject *)&PyTree_6dLType);
PyModule_AddObject(m, (char*)"KDTree_7Double", (PyObject *)&PyTree_7dLType);
PyModule_AddObject(m, (char*)"KDTree_8Double", (PyObject *)&PyTree_8dLType);
PyModule_AddObject(m, (char*)"KDTree_9Double", (PyObject *)&PyTree_9dLType);
PyModule_AddObject(m, (char*)"KDTree_10Double", (PyObject *)&PyTree_10dLType);
}


