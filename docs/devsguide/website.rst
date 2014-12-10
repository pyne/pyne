===========================
Building the PyNE website
===========================

-------------
Dependencies
-------------

Building the website/documentation requires the following:

   #. `Sphinx <http://sphinx-doc.org/>`_
   #. `sphinxcontrib-bibtex <https://pypi.python.org/pypi/sphinxcontrib-bibtex/>`_
   #. `PrettyTable <https://code.google.com/p/prettytable/>`_
   #. `Cloud Sphinx Theme <https://pythonhosted.org/cloud_sptheme/cloud_theme.html>`_

-----------------------------------
Procedure for modifying the website
-----------------------------------

The PyNE website source files are located in the ``docs`` directory. A developer first
makes necessary changes, then rebuild the website locally by executing the command::

    make html

This will generate html files for the website in the ``_build/html/`` folder.
The developer may view the local changes by opening these files with their 
favorite browser, e.g.::

    google-chrome _build/html/index.html

Once the developer is satisfied with the changes, the changes should be
commited and pull-requested per usual. Once the pull request is accepted, the
developer can push their local changes directly to the website by::

    make push-root

