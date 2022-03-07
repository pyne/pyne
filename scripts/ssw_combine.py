"""This example uses the SurfSrc class in pyne.mcnp to combine multiple files
together.

Script includes method combine_multiple_ss_files to combine MCNP surface source
entries with a combined header. The SurfSrc.__cmp__() method is reimplemented
to compare only the parts of the ssr files that matter for this operation.
"""

from itertools import izip

from numpy import copysign

from pyne.mcnp import SurfSrc


def combine_multiple_ss_files(newssrname, ssrnames):
    """Method reads headers from ssr1name and ssr2name binary files
    and checks if the headers are 'similar.'  If not, it returns False.

    Parameters
    ----------
    newssrname : str
        Path to new ssr file.
    ssrnames : list of str
        List of paths to ssr files being combined.

    Returns
    -------
    bool
        True if successfully created the new ssr file.

    """
    ssrfiles = [SurfSrc(ssrname, "rb") for ssrname in ssrnames]

    # read first surface source file's header
    ssrfiles[0].read_header()
    # then compare it with the headers from other files.
    for cnt, ssrfile in enumerate(ssrfiles[1:], 1):
        ssrfile.read_header()
        # we quit if there is a mismatch.
        if _compare_compatible(ssrfiles[0], ssrfile) != 0:
            print(
                "Headers do not match for all files.\nFile #{0} does not "
                "match 1st file.".format(ssrnames[cnt])
            )
            return False

    # calculate list of offsets for offsetting each track's nps value.
    trackoffsets = [
        sum(ssrfile.np1 for ssrfile in ssrfiles[:x]) for x in xrange(len(ssrfiles))
    ]

    ######################
    # Create new ssr file's header primarily from first ssr file's header
    newssr = SurfSrc(newssrname, "wb")

    # header
    newssr.kod = ssrfiles[0].kod
    newssr.ver = ssrfiles[0].ver
    newssr.loddat = ssrfiles[0].loddat
    newssr.idtm = ssrfiles[0].idtm
    newssr.probid = ssrfiles[0].probid
    newssr.aid = ssrfiles[0].aid
    newssr.knod = ssrfiles[0].knod

    # table 1
    newssr.np1 = sum(x.orignp1 for x in ssrfiles)  # note orignp1 vs np1
    newssr.nrss = sum(x.nrss for x in ssrfiles)
    newssr.ncrd = ssrfiles[0].ncrd
    newssr.njsw = ssrfiles[0].njsw
    newssr.niss = ssrfiles[0].niss
    newssr.table1extra = ssrfiles[0].table1extra

    # table 2
    newssr.niwr = ssrfiles[0].niwr
    newssr.mipts = ssrfiles[0].mipts
    newssr.kjaq = ssrfiles[0].kjaq
    newssr.table2extra = ssrfiles[0].table2extra

    # surfaces list
    newssr.surflist = ssrfiles[0].surflist

    # whatever the range from njsw to njsw+niwr is...
    for j in range(ssrfiles[0].njsw, ssrfiles[0].njsw + ssrfiles[0].niwr):
        print("unsupported entries; needs mcnp.py additions")

    # copy summary info/table
    newssr.summary_table = ssrfiles[0].summary_table
    newssr.summary_extra = ssrfiles[0].summary_extra

    # put the above header contents into the new ssr file
    newssr.write_header()

    #####################
    # Write each file's particle tracks in order of their listing.
    # The nps value is offset for each track by the entries in trackoffsets.
    ssrfiles[0].read_tracklist()
    newssr.tracklist = ssrfiles[0].tracklist

    for ssrfile, trackoffset in izip(ssrfiles[1:], trackoffsets[1:]):
        ssrfile.read_tracklist()
        for cnt, trackdata in enumerate(ssrfile.tracklist):
            trackdata.nps += copysign(trackoffset, trackdata.nps)
            ssrfile.tracklist[cnt] = trackdata
        newssr.tracklist.extend(ssrfile.tracklist)

    newssr.write_tracklist()

    print("Finished writing to new surface source file '{0}'" "".format(newssrname))
    newssr.close()

    return True


def _compare_compatible(first, other):
    """Determine if two ssr files are compatible for combining together.

    Parameters
    ----------
    first, other : SurfSrc objects
        Two ssr files to be compared.

    Returns
    -------
    int
        0 if all comparisons are valid, else the result of cmp() which was
        determined to be invald.
    """

    # Same code version/build used
    if other.kod != first.kod:  # code name
        return cmp(other.kod, first.kod)
    if other.ver != first.ver:  # major code version
        return cmp(other.ver, first.ver)
    if other.loddat != first.loddat:  # code version date
        return cmp(other.loddat, first.loddat)

    # Particle type
    if other.mipts != first.mipts:  # par type
        return cmp(other.mipts, first.mipts)

    # Geometry features
    if other.njsw != first.njsw:  # num of surfaces
        return cmp(other.njsw, first.njsw)
    if other.niwr != first.niwr:  # cells
        return cmp(other.niwr, first.niwr)
    if other.kjaq != first.kjaq:  # macrobody facet flag
        return cmp(other.kjaq, first.kjaq)
    for surf in range(len(first.surflist)):
        if other.surflist[surf].id != first.surflist[surf].id:
            # ID doesn't match
            return cmp(other.surflist[surf].id, first.surflist[surf].id)
        if other.surflist[surf].facet_id != first.surflist[surf].facet_id:
            # facet_id doesn't match
            return cmp(other.surflist[surf].facet_id, first.surflist[surf].facet_id)
        if other.surflist[surf].type != first.surflist[surf].type:
            # type doesn't match
            return cmp(other.surflist[surf].type, first.surflist[surf].type)
        if other.surflist[surf].num_params != first.surflist[surf].num_params:
            # num_params ddoesn't match
            return cmp(other.surflist[surf].num_params, first.surflist[surf].num_params)
        if other.surflist[surf].surf_params != first.surflist[surf].surf_params:
            # surf_params doesn't match
            return cmp(
                other.surflist[surf].surf_params, first.surflist[surf].surf_params
            )
    return 0


if __name__ == "__main__":
    import sys

    print(sys.argv)
    combine_multiple_ss_files("newssr", sys.argv[1:])
