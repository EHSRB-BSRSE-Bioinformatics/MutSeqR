# This file tells R CMD check about "global" variables that are used
# in the bundled PROAST code. This is necessary to avoid NOTES during
# the check process, as these PROAST variables are not used by MutSeqR

utils::globalVariables(
  names = c(
    # PROAST internal functions
    "f.adjust.saved", "f.cat", "f.ced.cat", "f.CI.sel.ord", "f.delete.gw",
    "f.dtype6.mn", "f.expect.cat", "f.firstcheck", "f.nested.con",
    "f.plot.cxt", "f.plot.sep", "f.press.key.to.continue", "f.quick.cat",
    "f.select.m46.con", "f.show.cat", "f.start.cat", "f.use.comma.for.decimal",

    # PROAST internal state variables
    ".gw.size", ".Pr.last", ".ypos"
  )
)