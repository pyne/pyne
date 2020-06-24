**Added:** in serpent.py, function VOL_serp2_fix() to correct for
*_VOL variable being an array. as seen in serpent 2

**Changed:** serpent.py function parse_dep.  Catches ValueError that
occurs when parse_dep attempts to make_mats with a serpent 2 *_dep.m file
due to the *_VOL variable not being a float

**Deprecated:** None

**Removed:** None

**Fixed:** issue #1224

**Security:** None
