.. _usersguide_jsoncpp:

=======================================
JsonCpp: Python bindings for JsonCpp
=======================================
The `JsonCpp project`_ provides an excellent in-memory JSON data structure as well 
as string writers and parsers.  Here are Python bindings for JsonCpp using Cython.
This module is mirror at the `PyJsonCpp project`_.

.. currentmodule: pyne.jsoncpp

.. _JsonCpp project: http://jsoncpp.sourceforge.net/

.. _PyJsonCpp project: https://github.com/scopatz/pyjsoncpp

-------------
Usage Example
-------------
The Value object is the main interface for 
Values may be converted to and from regular Python types.  These have the 
normal behavior for their type.

.. code-block:: python

    >>> from pyne.jsoncpp import Value, Reader, FastWriter, StyledWriter

    >>> v = Value({'name': 'Terry Jones', 'age': 42.0})
    >>> v['quest'] = "To find the grail."
    >>> v.keys()
    ['age', 'name', 'quest']
    >>> v['name']
    'Terry Jones'

    >>> v = Value([1, 2, 5, 3])
    >>> v[1:-1]
    [2, 5]
    >>> v[1:-1] = [42, 65]
    >>> v
    [1, 42, 65, 3]

    >>> v = Value("No one expects the Spanish Inquisition!!!!")
    >>> len(v)
    42

The Python Value class provides a view into the underlying C++ class.
This allows you to create several views into the same data.  For example, 
start with the following nested dictionary:

.. code-block:: python

    # make a nested dict and a view into a top-level item
    >>> v = Value({'a': {'b': 14}})
    >>> view_a = v['a']
    >>> view_a 
    {"b":14}

    # add an item to the view
    >>> view_a['c'] = 16
    >>> view_a 
    {"b":14,"c":16}

    # this item is present in the original value
    >>> v
    {"a":{"b":14,"c":16}}

Furthermore, there is a Reader class for converting JSON strings or files into 
Value instances.  There are also two writer classes, FastWriter and StyledWriter, 
for converting Value instances into compact and human-readable strings respectively.
For example:

.. code-block:: python

    >>> v = Value({'hello': 1})
    >>> fw = FastWriter()
    >>> fw.write(v)
    '{"hello":1}\n'

    >>> sw = StyledWriter()
    >>> print sw.write(v)
    {
       "hello" : 1
    }

    >>> r = Reader()
    >>> new_v = r.parse('{"hello":1}\n')
    >>> isinstance(new_v, Value)
    True
