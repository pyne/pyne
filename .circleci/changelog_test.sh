#!/bin/sh

# default main repo setup
master_repo="https://github.com/pyne/pyne.git"
default_branch="develop"

echo "*********"
echo "*********"
echo "*********"
echo "*********"
echo "MY WAR"
echo "*********"
echo "*********"
echo $GITHUB_PR_BASE_BRANCH
echo "*********"
echo "*********"
echo "*********"
echo "*********"


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
