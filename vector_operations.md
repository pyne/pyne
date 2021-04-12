Please follow these guidelines when making a Pull Request.
This template was adapted from [here](https://github.com/stevemao/github-issue-templates/blob/master/questions-answers/PULL_REQUEST_TEMPLATE.md) and [here](https://github.com/stevemao/github-issue-templates/blob/master/conversational/PULL_REQUEST_TEMPLATE.md).

## Description
This PR introduces the capability to perform arithmetic operations on vector-formatted tag data. There are also tests added to ensure functionality.

## Motivation and Context
Previously, the operations on tag data could only be performed on scalars, so this PR will enable us to perform these on any type of tag data.

## Changes
Feature: tag data of any format can be performed on

## Behavior
The current functionality of the _do_op() method within mesh.py is that it works on scalar data. This will change to any type of data, vector or scalar.

## Other Information
Added tests to the test_mesh.py file. 

## Changelog file
All pull requests are required to update the [CHANGELOG](https://github.com/pyne/pyne/blob/develop/CHANGELOG.rst) file with the PR.  Your update can take different forms including:

* creating a new entry describing your change, including a reference to this pull request number
* adding a reference to your pull request number to an exist entry if that entry can include your changes
