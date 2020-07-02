PyNE-CI Documentation
=====================


PyNE Circle-CI use CIrcle_CI version 2.1 and is composed of 4 main sections:
`executors`, `commands`, `jobs` and `workflow`.

Executors
---------
The `executors` define aliases to the environment to run the different jobs. In our case these are
different docker images containing different combinations of the optional
dependencies of PyNE (`MOAB`, `pyMOAB`, `DAGMC`) in one of the two supported
Python versions (Python 2.7 or Python 3.6).


Commands
---------
The `commands` section provides the definitions of macro commands with arguments to be used
in the different jobs.

`news_check`: checks the presence of a new news file in the news folder


`save_container`: saves the status of a container to be reloaded in an
other job.

    - `arguments`: 
        - `build` (string): build configuration parameters, used to identify a
          build configuration


`pull_container`: pulls a previously saved (using the `save_container` command) status of container.

    - `arguments`:
        - `build` (string): build configuration parameters (has to match the build
          argument used to save the container)


`checkout_build`: checks out PyNE branch to be tested, builds it and saves it
using the `save_container` command.
    
    - `arguments`:
        - `build` (string): used to identify a configuration
        - `flags` (string): flags to be use when building PyNE


`run_test`: pulls a previous container using `pull_container` command, then runs the PyNE
nosetests using the provided `flag`.
    
    - `arguments`:
        - `build` (string): used to identify a configuration
        - `flags` (string): flags to be use when running the PyNE nosetests


`website_build_push`: pull the Python 2.7 build with all the PyNE optional
depedencies `python2_dagmc_pymoab` saved image (using `pull_container` command),
build the website, and push it to the `pyne.github.com` repo.
    
    - `arguments`:
        - `push_option` (string: `test` or `root`): 
            - `test` option will push the newly built website to `website_preview` 
            branch of the repo allowing reviews. 
            - `root` option will push the website on the `master` branch of the repo, 
            deploying a new version of the website in `pyne.io`


Jobs
----
The `jobs` section defines all the different jobs used in the different
workflows. For each built configuration two jobs have been defined, one to build
PyNE, and one to run the `nosetests`. In addition to the build and test
jobs, two additional jobs have been created to build and push the website, one
to the `website_test` branch of the `pyne.github.com` repo, and one to the
`master` branch.


Workflow
--------
The `workflow` section defines independent workflow triggers the different
jobs (with triggers).
