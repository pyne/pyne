*********************
New Developer's Guide
*********************
=================
Setting up Github
=================
Before you start using GitHub, you have to make Git available on your computer.
Even if it’s already installed, it’s probably a good idea to update to the 
latest version. You can either install it as a package or via another installer, 
or download the source code and compile it yourself.

For download of Git visit `here <http://git-scm.com/downloads>`__

-----------------
Installing on Mac
-----------------
On Mavericks (10.9) or above you can do this simply by trying to run git from 
the Terminal the very first time. If you don’t have it installed already, it 
will prompt you to install it. You might be prompted to install command line 
developer tools if Xcode isn’t installed on your computer. If you want a more up 
to date version, you can also install it via a binary installer. An OSX Git 
installer is maintained and available for download at the `Git website 
<http://git-scm.com/download/mac.>`__.

---------------------
Installing on Windows
---------------------
There are also a few ways to install Git on Windows. The most official build is 
available for download on the Git `website <http://git-scm.com/download/win>`__ 
and the download will start automatically. Note that this is a project called 
Git for Windows (also called msysGit), which is separate from Git itself; for 
more information on it, go to http://msysgit.github.io/.

-------------------
Installing on Linux
-------------------
If you want to install Git on Linux via a binary installer, you can generally do 
so through the basic package-management tool that comes with your distribution. 
If you’re on Fedora for example, you can use yum:

$ sudo yum install git

If you’re on a Debian-based distribution like Ubuntu, try apt-get:

$ sudo apt-get install git

For more options, there are instructions for installing on several different 
Unix flavors on the Git website, at http://git-scm.com/download/linux.

==================================
Connecting Github account to Git
==================================
More detailed steps go here 
‘Github<https://help.github.com/categories/bootcamp/>’

#. In your computer's terminal, tell Git your name so your commits will be 
properly labeled. Type everything after the $ here:

$ git config --global user.name "YOUR NAME"

#. Tell Git the email address that will be associated with your Git commits. The 
email you specify should be the same one found in your email settings on Github. 
To keep your email address hidden click 'here 
<https://help.github.com/articles/keeping-your-email-address-private/>'_ 

$ git config --global user.email "YOUR EMAIL ADDRESS"

==================================
Installing from Source
==================================
NOTE: Might have to install pip, home-brew, and other programs to correctly 
install.
      Also, installing PyNE isn’t the same as cloning the repository.

PyNE's `website <http://pyne.io/install/index.html>`__ will lead you through the 
installation. 

To contribute to the project, 

#. Fork PyNE’s `Github repository <https://github.com/pyne/pyne>`__
 (Pull requests can be done on Github by comparing original to your forked 
repository). 

#. On Github, copy the HTTPS(recommended) link on your the page of your forked 
version that is located on the right side of the page. 
 
#. To keep things organized, it would be best to make a directory somewhere that 
is going to keep these files (Use “cd” to navigate to desired location and use 
“mkdir DESIRED NAME OF FOLDER” to make the directory and then go inside it). 
Make that directory your working directory and then enter ::

	$ git clone https://github.com/pyne/pyne.git

#. The PyNE files are now ready for your manipulations. Everything seen on the 
PyNE website is given to you. You can easily contribute by editing the contents 
of the folders, submitting these changes to your repository, and making a pull 
request to PyNE through Github’s website (click on the Pull Requests tab on the 
right side of the GitHub page and then submit a New Pull Request).

==================================
Signing up for list hosts 
==================================
Writing effective code isn’t easy. Thankfully, the PyNE developers can always be 
contacted on the listhost at pyne-dev@groups.google.com. Another way to get help 
is going to https://groups.google.com/forum/#!forum/pyne-users and joining the 
group to post. 

================
Getting Practice 
================
Novices to open-source projects can get still beneficially contribute to PyNE.  
To do so, go to PyNE’s ‘GitHub page <https://github.com/pyne/pyne>’ and, on the 
right hand side of the page, click on Issues. Once on this page, click on the 
“low hanging pinoli” label to display issues beginners can solve.

.. figure:: lhp.png :width: 700px
    :align: center
    :height: 340px

Also, if you were wondering, “low hanging pinoli” is a pine pun for low hanging 
pine fruit to call newbies. 

==================================
Adding and Updating Documentation 
==================================
To contribute, you can edit the text file in any program that allows you to edit 
text(Vim,textedit, Nano, etc) and doesn’t invisibly add characters to the 
file(like Word). The only important part is to write the file in a manner that’s 
considered reStructuredText (check out http://sphinx-doc.org/rest.html). Then, 
Sphinx will do everything else under the hood as described `here 
<http://pyne.io/devsguide/website.html>`__. Finally, commit these changes to 
your forked version and submit a pull request (through GitHub or the command 
line). 
