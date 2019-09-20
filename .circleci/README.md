PyNE-CI Documentation
=====================


PyNE Circle-CI use CIrcle_CI version 2.1 and is composed in 4 main sections:
`executors`, `commands`, `jobs` and `workflow`.

Executors
---------
The `executors` defines aliases to the environment to run the differents jobs, in our case
different docker images containing differents combinaisons of the optionnal
dependencies of PyNE (`MOAB`, `pyMOAB`, `DAGMC`) in one of the two supported
PyTHON flavor (Python 2.7 or Python 3.6).


Commands
---------
The `commands` section allows to define macro commands with arguments to be used
in the differents jobs.
`news_check`: check the presence of a new news file in the news folder

`save_container`: allows to save the status of a container to be reloaded in an
other job.
    `arguments`: 
        - `build` (string): build configuration parameters, used to identify a
          build configuration

`pull_container`: pulls a previously saved (using the `save_container` command) status of container.
    `arguments`:
        - `build` (string): build configuration parameters (has to match the build
          argument used to save the container)

`checkout_build`: checkout PyNE branch to be tested, build its and saves it
usins the `save_container` command.
    `arguments`:
        - `build` (string): used to identify a configuration
        - `flags` (string): flags to be use when building PyNE

`run_test`: pull a previous using `pull_container` command, then run the PyNE
nosetests using the provided `flag`.
    `arguments`:
        - `build` (string): used to identify a configuration
        - `flags` (string): flags to be use when running the PyNE nosetests

`website_build_push`: pull the Python 2.7 build with all the PyNE optionnal
depedencies `python2_dagmc_pymaob` saved image (using `pull_container` command)
build the website and push it to the `pyne.github.com` repo.
    `arguments`:
        `push_option` (string: `test` or `root`): push `test` will push the new
        built website to `website_preview` branch of the repo allowing previous
        allowing reviews. the push `root` version will push the website on the
        `master` branch of the repo, deploying a new version of the website in
        `pyne.io`

Jobs
----
The `jobs` section defines all the different jobs used in the different
workflows. For each built configuration two jobs have been defined, one to build
PyNE, and one to run the `nosetests`. In addition to the the build and test
jobs, two additional jobs have been created to build and push the website, one
to the `website_test` branch of the `pyne.github.com` repo, and one to the
`master` branch.


Workflow
--------
The `workflow` section defines independant workflow triggers the differents
jobs (with triggers).
