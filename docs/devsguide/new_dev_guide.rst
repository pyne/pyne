.. _devsguide_new_dev_guide:

*********************
New Developer's Guide
*********************

==================
Setting up Github:
==================
Before you start using GitHub, you have to make Git available on your computer.
Even if it’s already installed, it’s probably a good idea to update to the
latest version. You can either install it as a package or via another installer
or download the source code and compile it yourself.

To download Git visit 
https://git-scm.com/book/en/v2/Getting-Started-Installing-Git


-----------------
Installing on Mac
-----------------
On Mavericks (10.9) or above you can do this simply by trying to run git from
the Terminal the very first time. If you don’t have it installed already, it
will prompt you to install it. You might be prompted to install command line
developer tools if Xcode isn’t installed on your computer. 

If you want a more up to date version, you can also install it via a binary 
installer. An OSX Git installer is maintained and available for download at 
the http://git-scm.com/download/mac. This method will require Homebrew, which
may not already be on your device.


To Install Homebrew
```````````````````
This method uses Ruby, but there are other methods you can use `here <https://docs.brew.sh/Installation>`_ .

``To Install Ruby:``

`Installing Ruby <https://www.ruby-lang.org/en/documentation/installation/>`_

``Resources:``

#. `Learning Ruby <http://rubylearning.com/satishtalim/tutorial.html>`_
#. `Official FAQ <https://www.ruby-lang.org/en/documentation/faq/>`_

Run the commands ::

	$ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
	$ brew doctor
	$ brew tap brewsci/science

---------------------
Installing on Windows
---------------------
There are also a few ways to install Git on Windows. The most official build is
available for download on the `Git website <http://git-scm.com/download/win>`__
and the download will start automatically. Note that this is a project called
Git for Windows (also called msysGit), which is separate from Git itself; for
more information on it, go to http://msysgit.github.io/.


-------------------
Installing on Linux
-------------------
If you want to install Git on Linux via a binary installer, you can generally do
so through the basic package-management tool that comes with your distribution.
If you’re on Fedora, for example, you can use yum: ::

$ sudo yum install git

If you’re on a Debian-based distribution like Ubuntu, try apt-get: ::

$ sudo apt-get install git

For more options, there are instructions for installing on several different
Unix flavors on the Git website, at http://git-scm.com/download/linux.

--------------------------------
Connecting Github account to Git
--------------------------------
For more detailed steps go `here 
<https://help.github.com/categories/bootcamp/>`__.

Overview:

1. In your computer's terminal, tell Git your name so your commits will be
properly labeled. ::
    $ git config --global user.name "YOUR NAME"
2. Tell Git the email address that will be associated with your Git commits. The
email you specify should be the same one found in your email settings on Github. ::
    $ git config --global user.email "YOUR EMAIL ADDRESS"
3. To contribute to the project, Fork PyNE’s `Github repository 
<https://github.com/pyne/pyne>`__.

Learn how to keep your email address hidden `here
<https://help.github.com/articles/keeping-your-email-address-private/>`__.

================
Installing PyNE:
================

There are a couple of methods outlined for downloading PyNE on this site
under `Installation <https://pyne.io/install/index.html>`__. For
developers, it is recommended that you install PyNE from the source.

Make sure that your device has all of the dependencies downloaded. You 
can check this by typing the name of the program into your command line 
followed by "-v" (e.g. $ python -v). If you are just starting as a 
developer, likely, you will not have all of the necessary components.

A particularly important detail to pay attention to is the version of
each dependency. For example, PyNE is currently supported in conda
install from 2.7.\* - 3.6.*, and attempting to install with a default
version of 3.7.\* will result in errors. Setting your version to be
compatible before you begin will help ensure you are successful.
For a comprehensive list of all the software that PyNE depends on, examine
:ref:`Dependencies` .


-----------------------
Creating an Environment
-----------------------

If you already have Python installed through Anaconda, but it is not
compatible with the current version of PyNE, follow these steps to
create an environment for the conda install method.

#. Open a new shell in your terminal.
#. Run the command ::

	$ conda create -n py36 python=3.6 anaconda
        
	*  To run Python 2.7, or some other version, enter ::

			$ conda create -n py27 python=2.7 anaconda
	* "py27" serves as the name of the environment and "python=2.7"
			tells anaconda which version of Python to use in that 
			environment.
#. To enter the environment run ::

	$ conda activate py36

#. To leave enter ::
	
	$ conda deactivate

5. Before continuing with the PyNE installation, ensure that the Python version is correct.



For additional support with creating and managing environments,
documentation can be found
`here <https://docs.conda.io/projects/conda/en/latest/user-guide/
tasks/manage-environments.html>`__.


==========================
Signing up for list hosts:
==========================
Everyone faces challenges sometimes when writing effective code. Thankfully, the PyNE developers can always be 
contacted on the list host at pyne-dev@groups.google.com. Another way to get 
help is by going to https://groups.google.com/forum/#!forum/pyne-users and 
joining the group to post. 


========================
Preparing to Contribute:
========================

From the documentation of git and GitHub, the skills of forking and
cloning are especially important to have mastered before beginning any
contribution to the repository.

Before forking this or any other repository, engage SSH keys. This will
make the process of cloning and making the repository a remote straight
forward. The process to enable SSH keys on your device is found in the
GitHub documentation
`here <https://help.github.com/en/github/authenticating-to-github/connecting-
to-github-with-ssh>`__.

To fork, clone, and then make the original repository a remote follow
these steps.

#. Go to `PyNE <https://github.com/pyne/pyne>`__.
#. Select Fork, and then your account.
#. In the command line, enter ::

	$ git clone git@github.com:USERNAME/pyne.git 
	Replace USERNAME with your account's name.
#. Then add the original repository as a remote to your account, first
   enter the clone on your device with the command "cd pyne".
#. Complete the process by entering ::
	
	$ git remote add upstream git@github.com:pyne/pyne.git


The PyNE files are now ready for your contributions. You can easily contribute 
by editing the contents of the folders, submitting these changes to your 
repository, and making a pull request to PyNE through Github’s website. 


=============
Contributing:
=============

Follow the `Developer's Guide <https://pyne.io/devsguide/index.html>`__
for contributions to this site, and PyNE itself.

----------------
Getting Practice
----------------
Novices to open-source projects can get still contribute to PyNE.  
To do so, go to PyNE’s `GitHub Page <https://github.com/pyne/pyne>`__ and, 
on the right-hand side of the page, select Issues. Once on this page, select 
the “low hanging pinoli” label to display more issues with the same tag.
Pinoli is the Italian word for the Pine Nut, and this marker is the 
first place New Developers should look to contribute.

---------------------
Making a Pull Request
---------------------
The pull request (PR) should contain a description of the the changes, appropriate labels 
(such as bug, docs, test), projects (such as Transition to PyMOAB), and a reference
to the issue that led to the PR, which you can do by inserting '#issue_number' to 
the description.

A specific description, accurate labels, and appropriate
project selections will help other contributors to discover and review your work more easily.
Adding a reference to the issue the pull request will allow people to see the issue 
inspired it alongside any conversation about the issue. Before you make the pull request (PR), 
make sure that you include a news file named after the purpose of your PR 
following the template located in the `news` directory.

---------------------------------
Adding and Updating Documentation 
---------------------------------
To contribute, you can edit the text file in any program that allows you to edit 
text (Vim, TextEdit, Nano, etc.) and doesn’t invisibly add characters to the 
file (like Word). The only important part is to write the file in a manner that’s 
considered reStructuredText (check out http://sphinx-doc.org/rest.html). Then, 
Sphinx will do everything else under the hood as described `here 
<http://pyne.io/devsguide/website.html>`__.  

Your contributions will be more robust 
if they follow the formatting of other documents in the PyNE repository. As such, 
before you create or update a file, it is a good idea to skim through other PyNE 
documentation to see how they are formatted. Finally, commit these changes to your 
forked version and submit a pull request. 

.. _Dependencies:

==================
PyNE Dependencies:
==================

Something to note before you begin the installation process is the version of Python 
that you are using. PyNE is currently supported from Python 2 through Python 3.6 
(NOT 3.7, which is coming in future updates).


------------
PyNE Library
------------

PyNE is what is known as a library in Python (as with many languages, 
libraries can be built and installed to add functionalities to the language), 
and it is designed to assist in Nuclear Engineering projects. 

**Python:**

Python is a popular computer language for building applications, 
data science, and visualization. 

``To Install:``

A comprehensive Python package is Anaconda 
which provides a free installation for MacOS, Windows, and Linux. Follow the instructions 
for your specific platform `here <https://docs.anaconda.com/anaconda/install/>`_ if you would 
like to use Anaconda. Installing through Anaconda will result in many of the packages and 
dependencies being satisfied.

If you want to download Python without the rest of Anaconda, use 
the links in the resources below. Keep in-mind that PyNE is currently 
not supported on Python 3.7 or later, but you can 
:ref:`create an environment<Creating an Environment>` following the instructions above
for Anaconda users.

``Resources:``

#. `Python Download (Not Anaconda) <https://www.python.org/downloads/>`_
#. `Python Wiki <https://wiki.python.org/moin/>`_
#. `Beginner’s Guide <https://wiki.python.org/moin/BeginnersGuide>`_


**Ruby:**

Ruby is a programming language with a wide array of functionality, and a 
syntax similar to that of Python.

``To Install:``

`Installing Ruby <https://www.ruby-lang.org/en/documentation/installation/>`_

``Resources:``

#. `Learning Ruby <http://rubylearning.com/satishtalim/tutorial.html>`_
#. `Official FAQ <https://www.ruby-lang.org/en/documentation/faq/>`_


**GCC:**

GCC stands for GNU Compiler Collection, and GNU is a Unix-Like operating system 
that allocates resources and communicates with hardware. GNU is often used with 
a kernel called Linux. The Compiler Collection is one element that allows for 
PyNE to work as a Python library, even though it is built with several languages.

``To Install:``

`Installing GCC <https://gcc.gnu.org/install/index.html>`_

Homebrew Install:

If you have homebrew installed on your device, you can install GCC by entering ::

	$ brew install gcc

``Resources:``

#. `Documentation <https://gcc.gnu.org/onlinedocs/gfortran/#toc-Compiler-Characteristics-1>`_
#. `GCC FAQ <https://gcc.gnu.org/wiki/FAQ>`_
#. `GCC Glossary <https://gcc.gnu.org/wiki/GCC_glossary>`_


**CMake:**

CMake is designed to build, test, and compile packages programs across platforms.

``To Install:``

To download CMake, you can visit `Downloads <https://cmake.org/download/>`_.

If you have homebrew installed on your device, you can install GCC by entering ::

	$ brew install gcc

``Resources:``

#. `Wiki <https://gitlab.kitware.com/cmake/community/-/wikis/home>`_
#. `Tutorial <https://cmake.org/cmake/help/latest/guide/tutorial/index.html>`_


**Numpy:**

Numpy is a library in Python that adds support to matricies and arrays, 
as well as mathematical functions acting on those arrays.

``To Install:``

For the latest version using pip, in the command line enter ::

	$ pip install numpy

For the latest version using Conda, in the command line enter ::

	$ conda install numpy


``Resources:``

#. `Tutorial <https://numpy.org/learn/>`_
#. `Documentation <https://numpy.org/doc/stable/>`_


**SciPy:**

SciPy is a Python library that adds functions for linear algebra, optimization, 
integration, and other scientific or engineering tasks.

``To Install:``

`Installation <https://www.scipy.org/install.html>`_

``Resources:``

#. `Getting Started <https://www.scipy.org/getting-started.html>`_
#. `Tutorial <https://docs.scipy.org/doc/scipy/reference/tutorial/index.html>`_
#. `Documentation <https://www.scipy.org/docs.html>`_


**Cython:**

Cython is a compiler that helps in making C or C++ extensions for python.

``To Install:``

For the latest version, in the command line enter ::

	$ pip install cython

If you have homebrew installed on your device, you can install Cython by entering ::

	$ brew install cython

``Resources:``

#. `Wiki <https://github.com/cython/cython/wiki>`_
#. `User's Guide <https://cython.readthedocs.io/en/latest/src/userguide/index.html>`_
#. `Cython <https://cython.org>`_


**HDF5:**

HDF5 (the Hierarchical Data Format version 5) is a format that supports large, 
complex data in a file directory like structure similar to how you might with your computer.

``To Install:``

To install from source code, 
follow the instructions `here <https://www.hdfgroup.org/downloads/hdf5/source-code/>`_ .

For the latest version, in the command line enter ::

	$ brew install hdf5

``Resources:``

#. `Examples <https://portal.hdfgroup.org/display/HDF5/HDF5+Examples>`_
#. `Learning HDF5 <https://portal.hdfgroup.org/display/HDF5/Learning+HDF5>`_
#. `Known Problems <https://portal.hdfgroup.org/display/support/HDF5%201.12.0#knownprob>`_


**PyTables:**

PyTables is a package for managing large hierarchical datasets.

``To Install:``

For a variety of installation instructions, 
follow the instructions `here <http://www.pytables.org/usersguide/installation.html>`_ .

For the latest version, in the command line enter ::

	$ brew install tables

``Resources:``

#. `FAQ <http://www.pytables.org/FAQ.html>`_
#. `Tutorial <http://www.pytables.org/usersguide/tutorials.html>`_
#. `Project Pointers <http://www.pytables.org/project_pointers.html>`_


**LAPACK:**

LAPACK (Liner Algebra Package) is a software library for numerical liner algebra.

``To Install:``

Other installation methods can be found `here <http://www.netlib.org/lapack/#_software>`_

For the latest version, in the command line enter ::

	$ brew install lapack

``Resources:``

#. `FAQ <http://www.netlib.org/lapack/faq.html>`_
#. `User's Guide <http://www.netlib.org/lapack/lug/>`_


**BLAS:**

BLAS (Basic Linear Algebra Subroutines) coordinates operations on vectors and matricies.

``To Install:``

Other installation methods can be found `here <http://www.netlib.org/blas/#_software>`_

For the latest version, in the command line enter ::

	$ brew install openblas

``Resources:``

#. `Documentation <http://www.netlib.org/blas/#_documentation>`_


**Numexpr:**

Numexpr is a fast numerical evaluation tool for numpy, ensuring that 
expressions operating on arrays are faster and take up less memory.

``To Install:``

For the latest version, in the command line enter ::

	$ pip install numexpr

``Resources:``

#. `PyPi Project Homepage <https://pypi.org/project/numexpr/>`_
#. `Github Repository <https://github.com/pydata/numexpr>`_


-------
Website
-------

**Sphinx:** 

A python based documentation generator that allows files to be written into HTML, LaTeX, 
ePub, Texinfo, pages, and plain text. Sphinx uses reStructuredText, which is a very 
straight-forward markup language.

``To Install:``

For the latest version, in the command line enter ::

	$ pip install sphinx

``Resources:``

#. `Sphinx <https://www.sphinx-doc.org/en/master/>`_
#. `Tutorial <http://matplotlib.sourceforge.net/sampledoc/>`_
#. `reStructuredText Cheat Sheet <https://docutils.sourceforge.io/docs/user/rst/cheatsheet.txt>`_


**Sphinxcontrib-bibtex:**

An extension allowing Sphinx to interact with BibTeX.

``To Install:``

For the latest version, in the command line enter ::

	$ pip install sphinxcontrib-bibtex

``Resources:``

#. `Documentation <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_ 
#. `Known Issues and Workarounds <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html#known-issues-and-workarounds>`_
#. `Example <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/quickstart.html#minimal-example>`_


**PrettyTable:**

PrettyTable is a python library that adds a lot of versatility to table creation.

``To Install:``

For the latest version, in the command line enter ::

	$ pip install prettytable

``Resources:``

#. `Tutorial <https://code.google.com/archive/p/prettytable/wikis/Tutorial.wiki>`_


**Cloud Sphinx:**

Cloud is a Sphinx theme that PyNE uses to generate its 
HTML documentation (like this site).

``To Install:``

For the latest version, in the command line enter ::

	$ pip install cloud_sptheme

``Resources:``

#. `Documentation <https://cloud-sptheme.readthedocs.io/en/latest/>`_
#. `Source <https://foss.heptapod.net/doc-utils/cloud_sptheme>`_


**Jupyter:**

If you have downloaded Python through Anaconda Jupyter requirements should 
be satasfied, but it never hurts to make sure. 

You can check there version by entering ::

	$jupyter —-version

``To Install:``

If you are going to use Python 2:

For the latest version, in the command line enter ::

	$ pip install jupyter

If you are going to use Python 3:

For the latest version, in the command line enter ::

	$ pip3 install jupyter

``Resources:``

#. `Additional Installation Information <https://jupyter.readthedocs.io/en/latest/install.html>`_
#. `Project Documentation <https://jupyter.readthedocs.io/en/latest/projects/doc-proj-categories.html#deployment>`_
