#!/bin/sh

# default main repo setup
PR_BASE_BRANCH=develop
echo "Testing changelog against $PR_BASE_BRANCH branch"

master_repo="https://github.com/pyne/pyne.git"
default_branch=$PR_BASE_BRANCH

# setup temp remote 
git_remote_name=ci_changelog_`git log --pretty=format:'%h' -n 1`
git remote add ${git_remote_name} ${master_repo}
git fetch ${git_remote_name}

# diff against temp remote
added_changelog_entry=$((`git diff ${git_remote_name}/${default_branch} -- CHANGELOG.rst |wc -l`))

# cleaning temp remote
git remote remove ${git_remote_name}


# analysing the diff and returning accordingly
if [ $added_changelog_entry -eq 0 ]; then
    exit 1
fi
