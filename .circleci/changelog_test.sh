#!/bin/sh

# default main repo setup
PR_NUMBER=$(echo "$CIRCLE_PULL_REQUEST" | sed "s/.*\/pull\///")
API_GITHUB="https://api.github.com/repos/$CIRCLE_PROJECT_USERNAME/$CIRCLE_PROJECT_REPONAME"
PR_REQUEST_URL="$API_GITHUB/pulls/$PR_NUMBER"
PR_RESPONSE=$(curl "$PR_REQUEST_URL")
PR_BASE_BRANCH=$(echo $PR_RESPONSE | tr '\r\n' ' ' | jq -e '.base.ref' | tr -d '"')
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
