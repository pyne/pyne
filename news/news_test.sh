#!/bin/sh
# folder to look for addition
folder=$1

# default main repo setup
master_repo="https://github.com/pyne/pyne.git"
default_branch="develop"


# setup temp remote 
git_remote_name=ci_news_`git log --pretty=format:'%h' -n 1`
git remote add ${git_remote_name} ${master_repo}
git fetch ${git_remote_name}

# diff against temp remote
added_news_file=$((`git diff ${git_remote_name}/${default_branch} --name-only $folder |wc -l`))

# cleaning temp remote
git remote remove ${git_remote_name}


# analysing the diff and returning accordingly
if [ $added_news_file -eq 0 ]; then
    exit 1
fi
