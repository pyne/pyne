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
import pyrax
import os

def push_rackspace(fname, cred_file='rs.cred'):
    pyrax.set_credential_file(cred_file)
    cf = pyrax.cloudfiles
    with open(fname, 'rb') as f:
        fdata = f.read()
    obj = cf.store_object("pyne-data", fname, fdata)


pyrax.set_setting("identity_type", "rackspace")
pyrax.set_setting('region', 'ORD')
pyrax.set_credential_file('rs.cred')
cf = pyrax.cloudfiles
print "list_containers:", cf.list_containers()
print "get_all_containers:", cf.get_all_containers()
push_rackspace('prebuilt_nuc_data.h5')
