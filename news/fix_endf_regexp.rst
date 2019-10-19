**Added:** None

**Changed:** endf.Library._read_headers() and regular expressions in endf.pyx
    Removed regexps: CONTENTS_R, SPACE66_R, NUMERICAL_DATA_R
    Added regexps:   SPACEINT11_R
    Added methods:   _isContentLine(parts)

**Deprecated:** None

**Removed:** None

**Fixed:** Misidentification of descriptive text in (MF,MT)=(1,451) as contents lines.

**Security:** None
