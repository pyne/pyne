**Added:** None

**Changed:** Regular expressions in endf.pyx
    Removed regexps: CONTENTS_R, SPACE66_R, NUMERICAL_DATA_R
    Added regexps:   ELESSFLOAT_R, SPACEINT11_R
    Added methods:   _isContentLine(parts), _isDataLine(parts)

**Deprecated:** None

**Removed:** None

**Fixed:** Misidentification of descriptive text in (MF,MT)=(1,451) as contents lines.

**Security:** None
