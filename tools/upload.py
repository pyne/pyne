"""
# This is a simple tool to update the prebuilt_nuc_data.h5

In order to use this script you need to create an rs.cred file from the api information on the google drive doc.

Prerequisites:

 * PyNE
 * pyrax

Update procedure:

 * Switch to the appropriate branch and install PyNE based on it
 * run the following command in this folder:: 

 nuc_data_make --clean=1 --fetch-prebuilt False --make-open-only True -o prebuilt_nuc_data.h5

 * then run::
 
 python upload.py
 
 * The new prebuilt_nuc_data.h5 should now be on rackspace it may take 12-24 hours for this to propagate to all the CDN nodes.
"""
from __future__ import print_function
import pyrax
import os


def push_rackspace(fname, cred_file="rs.cred"):
    pyrax.set_credential_file(cred_file)
    cf = pyrax.cloudfiles
    with open(fname, "rb") as f:
        fdata = f.read()
    cont = cf.get_container("pyne-data")
    obj = cf.store_object("pyne-data", fname, fdata)
    cont.purge_cdn_object(fname)


pyrax.set_setting("identity_type", "rackspace")
pyrax.set_setting("region", "ORD")
pyrax.set_credential_file("rs.cred")
cf = pyrax.cloudfiles
print("list_containers: {}".format(cf.list_containers()))
print("get_all_containers: {}".format(cf.get_all_containers()))
push_rackspace("prebuilt_nuc_data.h5")
