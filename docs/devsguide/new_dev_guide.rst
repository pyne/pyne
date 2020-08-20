.. _devsguide_new_dev_guide:

*********************
New Developer's Guide
*********************

==================
Setting up Github:
==================

Before you start using GitHub, you have to install Git on your computer.
Even if it is on your device, it is a good idea to update to the latest version.
You can either install it as a package, via another installer or download the
source code and compile it yourself.

Download Git `here
<https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`__.

------------------------
Connecting GitHub to Git
------------------------
To connect Git on your device to GitHub, follow the instructions `here
<https://docs.github.com/en/github/getting-started-with-github/set-up-git#setting-up-git>`__.

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
followed by :code:`--version` (e.g., :code:`$ python --version`). If you are starting as a
developer, likely, you will not have all of the necessary components.


.. _creating_an_environment:

-----------------------
Creating an Environment
-----------------------

If you already have Python installed through Anaconda, but it is not
compatible with the current version of PyNE, follow these steps to
create an environment for the conda install method.

Documentation on creating and managing an environment can be found
`here <https://docs.conda.io/projects/conda/en/latest/user-guide/
tasks/manage-environments.html>`__.


==========================
Signing up for list hosts:
==========================

Everyone faces challenges sometimes when writing code. Thankfully, you can always
contact the PyNE Developers at pyne-dev@groups.google.com. Another way to
communicate with PyNE developers is to join the `Google Group
<https://groups.google.com/forum/#!forum/pyne-users>`__.


========================
Preparing to Contribute:
========================

The skills of forking and cloning are especially important to have mastered before
beginning any contribution to the repository.

Before forking this or any other repository, engage SSH keys. The process of
creating SSH keys on your device is detailed
`here <https://help.github.com/en/github/authenticating-to-github/connecting-
to-github-with-ssh>`__.

To fork, clone, and then make the original repository a remote follow
these steps.

#. Go to `PyNE <https://github.com/pyne/pyne>`__.
#. Select Fork, and then your account.
#. In the command line, enter ::

	$ git clone git@github.com:USERNAME/pyne.git 
	
        * Replace USERNAME with your account's name.
#. Then add the original repository as a remote to your account, first
   enter the clone on your device with the command "cd pyne".
#. Complete the process by entering ::
	
	$ git remote add upstream git@github.com:pyne/pyne.git

You can now easily contribute by editing the contents of the folders, pushing
these changes to your fork, and then making a pull request to the PyNE Github.


=============
Contributing:
=============

Follow the `Developer's Guide <https://pyne.io/devsguide/index.html>`__
for contributions to this site, and PyNE itself.

*Imposter syndrome disclaimer*: We want your help. No really, we do.

There might be a little voice inside that tells you you're not ready; 
that you need to do one more tutorial, or learn another framework, or 
write a few more of your own projects before you can contribute here.

We want to assure you, that's not the case.

This project has some clear Contribution Guidelines and expectations 
that you can read about in the Developer's Guide and below.

Below, you will find an outline of the process that you'll need to 
follow to get a pull request merged. By making expectations and process 
explicit, we hope it will make it easier for you to contribute.

And you don't just have to write code. You can help out by writing 
documentation, tests, or even by giving feedback about this work. 
(And yes, that includes giving feedback about the contribution guidelines.)

Thank you for contributing!

This disclaimer was originally written by `Adrienne Lowe 
<https://github.com/adriennefriend/imposter-syndrome-disclaimer/blob/master/README.md>`_ 
for a `PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted for
PyNE's Guide for New Developers.

----------------
Getting Practice
----------------
Novices to open-source projects can get still contribute to PyNE.  
To do so, go to PyNE’s `GitHub Page <https://github.com/pyne/pyne/issues>`__. Once
on this page, select the “low hanging pinoli” label to display more issues with the
same tag. Pinoli is the Italian word for the Pine Nut, and this marker is the
first place New Developers should look to contribute.

---------------------
Making a Pull Request
---------------------
Before you make the pull request (PR), make sure that you include a change to the 
CHANGELOG.rst file under the appropriate heading. Your contribution can be as simple 
as a line of text, but it should incorporate the PR's number in parenthesis at the end. 
If there is a line already in the changelog that describes a similar action (e.g., 
adding a publication) to your own, you can add the number of your PR after the existing line.

The PR should contain a description of the changes, appropriate labels, projects,
and a reference to the issue that led to it, which you can do by inserting
‘#issue_number’ to the description (i.e., adding #161 to your PR description).
These elements will help other contributors discover and review your work. Adding a
reference to the issue the pull request will allow people to see the issue inspired
it alongside any conversation about the issue.

----------------------
Changing Documentation
----------------------
To contribute, you can edit the text file in any program that allows you to edit
text (Vim, TextEdit, Nano, etc.) and does not invisibly add characters to the file
(like Word).

Your contributions will be more robust if they follow the formatting of other
documents in the PyNE repository. As such, before you create or update a file, it
is a good idea to skim through other PyNE documentation to see how they are
formatted. Finally, commit these changes to your forked version and submit a pull
request.
