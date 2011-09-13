.. currentmodule:: pyne.stlconverters

.. _usersguide_stlconverters:

=============================
C++ Standard Library Wrappers
=============================
Because they are useful, PyNE exposses a Python interface to many often used C++ standard libray
containers.  These are primarly used to efficiently deliver data from the low- to high-level
without excessive copying.

These wrapper classes allow you to initialize and/or modify the C++ containters after which they 
are named.  Thus the wrappers give you a significant amount of control over the memory.  You can 
choose whether or not you want to make a new instance of the underlying class.  (When using this 
module from Python it is *highly* recommended that you always make a new object.  From Cython you 
may not want to because you already have a pointer to the object you want to wrap.)  Additionally, 
you may also choose whether, on the dereferncing of the wrapper object, if you also want to 
deallocate the pointer.  (You may not want to deallocate if the pointer is shared by multiple wrappers.)  

These wrappers exist in a separate Cython module. Please feel free to compile and link against them.
Or if you need a wrapper for a type that is not included, these wrappers make great templates
for your code as well!

--------------
Example of Use
--------------
From Python, the following represent common use cases::

    import pyne.stlconverters as conv

    # New integer set
    s = conv.SetInt()
    s.add(7)
    assert (7 in s)
    assert (11 not in s)

    # New string set
    s = conv.SetStr(["Aha", "Take", "Me", "On"])
    assert ("Aha" in s)
    assert ("Captain Hammer" not in s)

    # Two new mapping from strings to ints
    m = conv.MapStrInt({'yes': 1, 'no': 0})
    assert (len(m) == 2)
    assert (m['no'] == 0)

    # Careful! We are only copying a view...
    n = conv.MapStrInt(m, False)
    assert (len(n) == 2)
    assert_equal(n['yes'] == 1)

    # ...so m & n point to the same underlying map!
    n['maybe'] = -1
    assert_equal(m['maybe'] == -1)

In Cython, the use case is a little different because we have access to pointers
on the C++ level.  Suppose we already have a map that exists, we simply want to wrap it
and expose the python bindings.

.. code-block:: cython

    cimport pyne.stlconverters as conv
    import pyne.stlconverters as conv

    # Existing map
    def conv._MapIntDouble i_exist_proxy = conv.MapIntDouble(False)
    i_exist_proxy.map_ptr = &i_exist_in_c
    i_exist_in_python = i_exist_proxy

    # Existing map & owned by someone else
    def conv._MapIntDouble owned_elsewhere_proxy = conv.MapIntDouble(False, False)
    owned_elsewhere_proxy.map_ptr = &owner.map_in_c
    owned_elsewhere_as_seen_in_python = owned_elsewhere_proxy

Further information on this wrapper module, please refer to the library reference :ref:`pyne_stlconverters`.
