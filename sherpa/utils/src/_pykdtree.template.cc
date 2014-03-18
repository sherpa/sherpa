
// %(STRUCTNAME)s
// %(TREENAME)s
// %(RECORDNAME)s
// %(POINTNAME)s
// %(COORDNAME)s


typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  %(TREENAME)s *tree;
} %(STRUCTNAME)s;


static PyObject* %(STRUCTNAME)s_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  %(STRUCTNAME)s *self;
  
  self = (%(STRUCTNAME)s*)type->tp_alloc(type, 0);
  if(self)
    self->tree = NULL;
  
  return (PyObject*)self;
}


static int %(STRUCTNAME)s_init(%(STRUCTNAME)s *self, PyObject *args, PyObject *kwds)
{
  self->tree = (%(TREENAME)s *) new %(TREENAME)s();
  return 0;
}


// New datatype dealloc method
static void %(STRUCTNAME)s_dealloc(%(STRUCTNAME)s* self) {
  if(self->tree)
    delete self->tree;
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* %(STRUCTNAME)s_repr(%(STRUCTNAME)s* self) {

  return PyString_FromString( (char*) "<%(CLASSNAME)s instance>" );
}

static PyObject* %(STRUCTNAME)s_str(%(STRUCTNAME)s* self) {

  std::ostringstream output;
  std::vector<%(RECORDNAME)s> *temp = 0;
  size_t size;
  size_t full = 6, half = 3;

  if(self && self->tree) {

    size = (int) self->tree->size();

    temp = (std::vector<%(RECORDNAME)s > *) self->tree->get_all();
    
    if (size > full) {
      
      std::vector<%(RECORDNAME)s >::const_iterator iter = temp->begin();
      for (size_t i = 0; i < half; i++, iter++)
	output << (*iter) << std::endl;
      
      output << "..." << std::endl;
      
      for (size_t i = size-half-1; i < size; i++) 
	output << temp->at(i) << std::endl;
      
    }

    else if( size <= full ) {
    
      std::vector<%(RECORDNAME)s >::const_iterator iter = temp->begin();
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



static PyObject * %(STRUCTNAME)s_add(%(STRUCTNAME)s* self, PyObject *args, PyObject *kwds) {

  %(RECORDNAME)s record;
  %(TREENAME)s *tree;

  if ( !PyArg_ParseTuple(args, "((%(DIMIDS)s)%(INDEXID)s)",
			 %(RECORDPOINTS)s
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (%(II)s dim %(DIMTYPE)s vector, %(INDEXTYPE)s value)");
    return NULL;
  }
  
  if(self && ((%(STRUCTNAME)s*)self)->tree) {

    tree = (%(TREENAME)s *)((%(STRUCTNAME)s*)self)->tree;
    tree->add(record);

  } else {
    PyErr_SetString(PyExc_TypeError,"Adding record failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}

static PyObject * %(STRUCTNAME)s_remove(%(STRUCTNAME)s* self, PyObject *args, PyObject *kwds) {

  %(RECORDNAME)s record;
  bool success;
  PyObject *result = Py_False;

  if ( !PyArg_ParseTuple(args, "((%(DIMIDS)s)%(INDEXID)s)",
			 %(RECORDPOINTS)s
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (%(II)s dim %(DIMTYPE)s vector, %(INDEXTYPE)s value)");
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

static PyObject * %(STRUCTNAME)s_size(%(STRUCTNAME)s* self) {

  int result;

  if(self && self->tree) {
    result = (int) self->tree->size();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing size() failed!");
    return NULL;
  }

  return Py_BuildValue((char*)"i", result);
}

static PyObject * %(STRUCTNAME)s_optimize(%(STRUCTNAME)s* self) {

  if(self && self->tree) {
    self->tree->optimize();
  } else {
    PyErr_SetString(PyExc_TypeError,"Accessing optimize() failed!");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject * %(STRUCTNAME)s_find_exact(%(STRUCTNAME)s* self, PyObject *args, PyObject *kwds) {

  %(RECORDNAME)s record;
  %(RECORDNAME)s *temp = NULL;
  PyObject *result = NULL;

  if ( !PyArg_ParseTuple(args, "((%(DIMIDS)s)%(INDEXID)s)",
			 %(RECORDPOINTS)s
			 &record.data) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (%(II)s dim %(DIMTYPE)s vector, %(INDEXTYPE)s value)");
    return NULL;
  }
  
  if(self && self->tree) {

    temp = (%(RECORDNAME)s *) self->tree->find_exact(record);

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
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(%(DIMIDS)s)",
						 %(TEMPPOINTS)s))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("%(INDEXID)s", temp->data))==-1) {
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



static PyObject * %(STRUCTNAME)s_find_nearest(%(STRUCTNAME)s* self, PyObject *args, PyObject *kwds) {

  %(POINTNAME)s point_t;
  %(COORDNAME)s *coord_t = NULL;
  %(RECORDNAME)s *temp = NULL;
  PyObject *result = NULL;
  %(TREENAME)s *tree;

  if ( !PyArg_ParseTuple(args, "(%(DIMIDS)s)",
			 %(POINTS1)s
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have %(II)s elements: %(II)s dim %(DIMTYPE)s vector");
    return NULL;
  }
  
  if(self && ((%(STRUCTNAME)s*)self)->tree) {

    tree = (%(TREENAME)s *)((%(STRUCTNAME)s*)self)->tree;

    coord_t = point_t;

    temp = (%(RECORDNAME)s *) tree->find_nearest(coord_t);

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
    
    if (PyTuple_SetItem(result, 0, Py_BuildValue("(%(DIMIDS)s)",
						 %(TEMPPOINTS)s))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(a) when setting element");
      
      Py_DECREF(result);
      delete temp;
      return NULL;
    }
    
    if (PyTuple_SetItem(result, 1, Py_BuildValue("%(INDEXID)s", temp->data))==-1) {
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

static PyObject * %(STRUCTNAME)s_count_within_range(%(STRUCTNAME)s* self, PyObject *args, PyObject *kwds) {

  %(POINTNAME)s point_t;
  %(COORDNAME)s *coord_t = NULL;
  size_t temp;
  PyObject *result = NULL;
  %(TREENAME)s *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(%(DIMIDS)s)d",
			 %(POINTS)s
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (%(II)s dim %(DIMTYPE)s vector, double value)");
    return NULL;
  }
  
  if(self && ((%(STRUCTNAME)s*)self)->tree) {

    tree = (%(TREENAME)s *)((%(STRUCTNAME)s*)self)->tree;

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



static PyObject * %(STRUCTNAME)s_find_within_range(%(STRUCTNAME)s* self, PyObject *args, PyObject *kwds) {

  %(POINTNAME)s point_t;
  %(COORDNAME)s *coord_t = NULL;
  std::vector<%(RECORDNAME)s> *temp = 0;
  %(TREENAME)s *tree;
  double val3;
  RANGE_T range;

  if ( !PyArg_ParseTuple(args, "(%(DIMIDS)s)d",
			 %(POINTS)s
			 &val3
			 ) ) {
    PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements: (%(II)s dim %(DIMTYPE)s vector, double value)");
    return NULL;
  }
  
  if(self && ((%(STRUCTNAME)s*)self)->tree) {

    tree = (%(TREENAME)s *)((%(STRUCTNAME)s*)self)->tree;

    coord_t = point_t;

    range = static_cast<RANGE_T >(val3);

    temp = (std::vector<%(RECORDNAME)s > *) tree->find_within_range(coord_t, range);

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

  std::vector<%(RECORDNAME)s >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(%(DIMIDS)s)%(INDEXID)s",
	%(ITERPOINTS)s, (*iter).data))==-1) {
      PyErr_SetString(PyErr_Occurred(),"(c) when setting element");
      
      Py_DECREF(result);
      return NULL;
    }
  }

  // delete the pointer after use!
  delete temp;
  
  return result;
}


static PyObject * %(STRUCTNAME)s_get_all(%(STRUCTNAME)s* self) {

  std::vector<%(RECORDNAME)s> *temp = 0;

  if(self && self->tree) {

    temp = (std::vector<%(RECORDNAME)s > *) self->tree->get_all();

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

  std::vector<%(RECORDNAME)s >::const_iterator iter = temp->begin();

  for (size_t i = 0; i < temp->size(); i++, iter++) {
    if (PyList_SetItem(result, i, Py_BuildValue("(%(DIMIDS)s)%(INDEXID)s",
	%(ITERPOINTS)s, (*iter).data))==-1) {
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

static PyMemberDef %(STRUCTNAME)s_members[] = {
  // {(char*)"name", T_OBJECT_EX, offsetof(%(STRUCTNAME)s, name), 0,
  //  (char*)"name"},

  // {(char*)"datatype", T_INT, offsetof(%(STRUCTNAME)s, datatype), 0,
  //  (char*)"datatype"},

  // {(char*)"unit", T_OBJECT_EX, offsetof(%(STRUCTNAME)s, unit), 0,
  //  (char*)"unit"},

  // {(char*)"desc", T_OBJECT_EX, offsetof(%(STRUCTNAME)s, desc), 0,
  //  (char*)"desc"},
  {NULL}  // Sentinel 
};

// Methods of class are inserted and 'named' here  

// Notice that the function signature type is mapped with a flag 
// METH_NOARGS, METH_VARARGS, or METH_VARARGS | METH_KEYWORDS 

static PyMethodDef %(STRUCTNAME)s_methods[] = {

  { (char*) "size",(PyCFunction)%(STRUCTNAME)s_size, METH_NOARGS, (char*) ""},

  { (char*) "optimize",(PyCFunction)%(STRUCTNAME)s_optimize, METH_NOARGS, (char*) ""},

  { (char*) "get_all",(PyCFunction)%(STRUCTNAME)s_get_all, METH_NOARGS, (char*) ""},

  { (char*) "add",(PyCFunction)%(STRUCTNAME)s_add, METH_VARARGS, (char*) ""},

  { (char*) "find_exact",(PyCFunction)%(STRUCTNAME)s_find_exact, METH_VARARGS, (char*) ""},

  { (char*) "find_nearest",(PyCFunction)%(STRUCTNAME)s_find_nearest, METH_VARARGS, (char*) ""},

  { (char*) "find_within_range",(PyCFunction)%(STRUCTNAME)s_find_within_range, METH_VARARGS, (char*) ""},

  { (char*) "count_within_range",(PyCFunction)%(STRUCTNAME)s_count_within_range, METH_VARARGS, (char*) ""},

  { (char*) "remove",(PyCFunction)%(STRUCTNAME)s_remove, METH_VARARGS, (char*) ""},

  {NULL, NULL, 0, NULL}  // Sentinel 
};


static PyTypeObject %(STRUCTNAME)sType = {
    PyObject_HEAD_INIT(NULL)
    0,                                // ob_size
    (char*)"%(CLASSNAME)s",           // tp_name, __name__
    sizeof(%(STRUCTNAME)s),           // tp_basicsize
    0,                                // tp_itemsize
    (destructor)%(STRUCTNAME)s_dealloc, // tp_dealloc
    0,                                // tp_print
    0,                                // tp_getattr, __getattr__
    0,                                // tp_setattr, __setattr__
    0,                                // tp_compare
    (reprfunc)%(STRUCTNAME)s_repr,    // tp_repr, __repr__
    0,                                // tp_as_number
    0,                                // tp_as_sequence
    0,                                // tp_as_mapping
    0,                                // tp_hash
    0,                                // tp_call, __call__
    (reprfunc)%(STRUCTNAME)s_str,     // tp_str, __str__
    0,                                // tp_getattro
    0,                                // tp_setattro
    0,                                // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags
    (char*)"%(STRUCTNAME)s objects",  // tp_doc, __doc__
    0,		                      // tp_traverse 
    0,		                      // tp_clear 
    0,		                      // tp_richcompare 
    0,		                      // tp_weaklistoffset 
    0,		                      // tp_iter, __iter__
    0,		                      // tp_iternext 
    %(STRUCTNAME)s_methods,           // tp_methods 
    %(STRUCTNAME)s_members,           // tp_members 
    0,                                // tp_getset 
    0,                                // tp_base 
    0,                                // tp_dict, __dict__
    0,                                // tp_descr_get 
    0,                                // tp_descr_set 
    0,                                // tp_dictoffset 
    (initproc)%(STRUCTNAME)s_init,    // tp_init, __init__
    0,                                // tp_alloc 
    %(STRUCTNAME)s_new,               // tp_new 
};

