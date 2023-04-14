"""JsonCpp tests"""
import os

try:
    import simplejson as json
except ImportError:
    import json
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


import glob


def test_notebooks():
    # add all notbooks with invalid JSON to this list
    failed = []

    # list all notebook files
    notebooks = glob.glob("../tutorial/*ipynb")
    for notebook in notebooks:
        try:
            with open(notebook, "r", encoding="utf-8") as myfile:
                try:
                    # try to parse the JSON
                    data = json.load(myfile)
                except:
                    failed.append(notebook)
        except TypeError:
            with open(notebook, "r") as myfile:
                try:
                    # try to parse the JSON
                    data = json.load(myfile)
                except:
                    failed.append(notebook)

    assert (
        len(failed) ==
        0,
        "Not all notebooks contained valid JSON. [%s] failed." % ", ".join(failed),
    )

