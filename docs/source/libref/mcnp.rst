.. _pyne_mcnp:

====================================================
MCNP Input and Output Interfaces -- :mod:`pyne.mcnp`
====================================================

.. currentmodule:: pyne.mcnp

.. automodule:: pyne.mcnp
  
There is a class for a variety of types of files that MCNP produces.  The
functionality of the module can be obtained by importing as such::

    from pyne import mcnp

*********************
Current functionality
*********************

************************
Functionality to discuss
************************

Anthony Scopatz was saying that he imagines this wrapper would operate through
dictionaries, somewhat along the idea of what json files are. The benefit is
mostly that the data is then not attached to a class. What I am interested in
is the object-oriented interface, which it seens can easily be a layer on top
of something like a json file. An NJOY project we found used json-type files to 

Perhaps the biggest objective with this code is creating something that is
persistent so that the input object is not created again every time the input
is to be modified. That is, an input object built from the interactive command
prompt could be loaded again without having to execute all the commands;
everything would be saved in the json file. Basically it seems that the
mcnpcard class would write to the json file instead. Anthony I want you to know
I think what you were saying now makes more sense to me and I think it is a
fantastic idea.

Error-checking: it would be nice to be able to do some of the error checking in
this class, especially in cases where the MCNPX errors are unclear. This may
not be relevant here, but for example the MCNPX error that results when the
line endings are not correct (LF versus CRLF) is not descriptive. Some of
MCNPX's error checking seems to be pretty good though.

Lots of repitition, and that's an issue I'm not clear how to solve: 
One of the things that is not great about what I've done is that for each card
there are two methods: one in mcnp.Inp and one in mcnpcard. The former has
information about the input as a whole, and the latter has the "low-level"
information about how the card is supposed to be written

Learning curve: What do we expect the user knows? When they use a temp keyword
on a cell card should the documentation tell them what free gas thermal
treatment is? Do we want people to have learned how to write cards before using
this wrapper, or be able to only know the wrapper interface?

**************
Inp Class
**************
.. autoclass:: Inp
   :members:
   :inherited-members:


.. automodule:: mcnpcard

