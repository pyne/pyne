**Added:** None

**Changed:** None

**Deprecated:** None

**Removed:** None

**Fixed:**

* Fixed issue where some gamma x-rays where throwing ``NotANuclide`` errors
  because the underlying nuclides were being read & recorded with negative ids.
  All nuclide ids are now ensured to be positive.

**Security:** None
