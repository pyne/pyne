Direct Accelerated Geometry Monte Carlo (DAGMC) Toolkit
==========================================================

The Direct Acclerated Geometry Monte Carlo (DAGMC) Toolkit is an
interface (`DAGMC Source
<http://trac.mcs.anl.gov/projects/ITAPS/browser/MOAB/trunk/tools/dagmc>`_)
to the `MOAB mesh database
<http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_ that provides the
methods necessary for ray tracing on a CAD-based geometric model.

This repository provides implementations of that interface for various
Monte Carlo radiation transport codes.  Depending on the distribution
limitations and software design of each code, there are many modes by
which we are delivering this capability.  There are 3 main
characteristics of a Monte Carlo physics package which define our
method of distribution.

Efforts are underway to make DAGMC available in the following physics
packages:
   * MCNP5: complete and in production use
   * Fluka: just beginning (12/2012)
   * Serpent: underway (12/2012)
   * OpenMC: planned for 2013
   * GEANT4: planned for 2013
   * MCNP6: planned for 2013

Geometry Interface
-------------------

In cases where the physics package has a cleanly defined geometry
interface (FLUKA), we are able to distribute a standalone collection of
methods that each user can compile and link with the physics package.

When the geometry interface is not cleanly defined, our modifications
include modifications to the original source code.  Therefore our
distribution mechanism depends the ability to integrate our
modifications into the main physics code development path.

Mainline Development Integration
----------------------------------

In cases where the authors of the physics package are willing to
integrate DAGMC as an option in their primary software distribution
(SHIFT), then the main distribution mechansim will be as part of that
software in its normal distribution channels.  This may either in be
in a shared software repository or in a regular release snapshot.

In cases where the primary authors prefer DAGMC to be distributed
separately as an enhancement for their software (MCNP5), the
distribution mechanism will be as a patch to their source code, since
we generally are not authorized to redistribute their code.

Building Documentation
-------------------------

(Note: Sphinx versions higher than 1.1.2 are required.)

A 2 branch system has been implemented to maintain a clean process of
rebuilding this site.

1. The `master` branch contains the restructured text documents and
Sphinx configuration used to build the site.  All direct editing of
files should be made in the `master` branch.

2. The `gh-pages` branch contains the processed and published web
content that is derived by Sphinx from the `master` branch.  These
files should not be editted directly.

Best practice workflow for contributing to site changes
--------------------------------------------------------

1. Checkout the `master` branch

   ```git checkout master```

2. Synchronize your branch with the repository (either `pull` or
`fetch` and `merge`)

     ```git pull upstream```

3. Create a branch to contain your change

     ```git checkout -b add_some_info```

4. Make your changes in this branch

5. Test you changes by using the `gh-preview` target

     ```make gh-preview```

   This will build a version of the site in the `gh-build` directory of
   your branch, `add_some_info`.  You can load it directly in a local
   browser.

6. Repeat steps 4-5 until satisfied.

7. Once satisfied with the master RST files, push your branch to the
   repo.  Be sure to synchronize with any possible changes to the
   `master` branch first.

     ```
     git fetch upstream
     git rebase upstream/master
     git push upstream add_some_info
     ```

8. Issue a pull request by going to your branch on the repo and
   clicking the "Pull Request" button.

Best practice for managing a pull request
------------------------------------------

1. Synchronize your repository with the remote repo

     ```git fetch upstream```

2. Checkout the `pull_request_branch`

     ```git checkout -b pull_request_branch upstream/pull_request_branch```

3. Test the changes by using the `gh-preview` target

    ```make gh-preview```

   This will build a version of the site in the `gh-build` directory in
   your branch, `pull_request_branch`.  You can load it directly in a
   local browser.

5. If satisfied, merge the `pull_request_branch` into the `master`
   branch.  Be sure to synchronize with the remote repo first.

     ```
     git checkout master
     git fetch upstream
     git rebase upstream/master
     git merge pull_request_branch
     ```

6. If there are no conflicts, push this to the repo

     ```git push upstream master```

7. Republish the pages with the `gh-publish` target.

     ```make gh-publish```

