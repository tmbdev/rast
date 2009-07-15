// -*- C++ -*-

%module rast
%include "typemaps.i"
%include "cstring.i"
%include "std_wstring.i"
%{
#include "rast.h"
#include "cedges.h"
%}

%inline %{
const int epsilon = 0;
%}

%exception {
    try {
        $action
    }
    catch(const char *s) {
        PyErr_SetString(PyExc_IndexError,s);
        return NULL;
    }
    catch(...) {
        PyErr_SetString(PyExc_IndexError,"unknown exception in iulib");
        return NULL;
    }
}

%{
int zero() { return 0; }
%}

%include "rast.h"
%include "cedges.h"
