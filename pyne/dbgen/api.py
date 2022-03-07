import tables as tb

build_dir = "build_nuc_data"

BASIC_FILTERS = tb.Filters(complevel=5, complib="zlib", shuffle=True, fletcher32=False)
