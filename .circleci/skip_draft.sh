#!/bin/bash

required_env_vars=(
    "GITHUB_TOKEN"
    "CIRCLE_PROJECT_USERNAME"
    "CIRCLE_PR_REPONAME"
    "CIRCLE_PR_NUMBER"
    "CIRCLE_TOKEN"
    "CIRCLE_BUILD_NUM"
)
for required_env_var in ${required_env_vars[@]}; do
    if [[ -z "${!required_env_var}" ]]; then
    printf "BAD ENVIRONMENT SETTINGS\n"
    printf "${required_env_var} not provided, but that doesn't mean we should skip CI.\n"
    exit 0
    fi
done
# Since we're piggybacking off of an existing OAuth var, tweak the var for our uses
token=$(printf "${GITHUB_TOKEN}" | cut -d':' -f1)
headers="Authorization: token $token"
api_endpoint="https://api.github.com/repos/${CIRCLE_PROJECT_USERNAME}/${CIRCLE_PR_REPONAME}/pulls/${CIRCLE_PR_NUMBER}"
# Fetch PR metadata from Github's API and parse fields from json
github_res=$(curl --silent --header "${headers}" "${api_endpoint}" | jq '{mergeable_state: .mergeable_state, title: .title}')
mergeable_state=$(printf "${github_res}" | jq '.mergeable_state')
title=$(printf "${github_res}" | jq '.title' | tr '[:upper:]' '[:lower:]')
echo "${title}"
if [[ "${title}" == "null" && "${mergeable_state}" == "null" ]]; then
    printf "Couldn't fetch info on PR, but that doesn't mean we should skip CI.\n"
    exit 0
fi

# Check for Draft
cancel_running_jobs=0
if [[ "${mergeable_state}" == "\"draft\"" ]]; then
    printf "PR is a draft, skipping CI!\n"
    cancel_running_jobs=1
fi
# Check for Title SKIP Token
for skip_token in '[skip ci]' '[ci skip]' '[wip]' '[WIP]'; do
    if [[ ${title} == *"${skip_token}"* ]]; then
    printf "Found \"${skip_token}\" in PR title, skipping CI!\n"
    cancel_running_jobs=1
    fi
done

# Cancel the Job if needed
if [[ "${cancel_running_jobs}" == 1 ]]; then
    printf "Attempting to cancel any running jobs"
    FORMAT_PARAM="-H \"Content-Type: application/json\""
    AUTH_PARAMS="--header  \"Circle-Token: ${CIRCLE_TOKEN}\""
    CMD="curl -X"
    CI_API_BASE_URL="https://circleci.com/api/v2"
    CI_PRJ_URL="${CI_API_BASE_URL}/project/gh/${CIRCLE_PROJECT_USERNAME}/${CIRCLE_PROJECT_REPONAME}"
    CI_JOB_URL="${CI_PRJ_URL}/job/${CIRCLE_BUILD_NUM}"
    # Get the Id of the Workflow to Cancel
    MSG=$(${CMD} GET "${CI_JOB_URL}" $FORMAT_PARAM -H "Circle-Token: ${CIRCLE_TOKEN}")
    WKF_ID=$(echo $MSG | tr '\r\n' ' ' | jq  ".latest_workflow.id" | tr -d '"')

    # Cancel the Workflow based on its Id
    printf "Canceling WORKFLOW: build_and_test"
    CI_WKF_CANCEL_URL="${CI_API_BASE_URL}/workflow/${WKF_ID}/cancel"
    curl -X POST "${CI_WKF_CANCEL_URL}" -H "Circle-Token: ${CIRCLE_TOKEN}"
else
    printf "No reason to skip CI, let's go!"
fi
exit 0
