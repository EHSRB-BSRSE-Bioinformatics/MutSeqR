# This file tells R CMD check about "global" variables that are used
# in the bundled PROAST code. This is necessary to avoid NOTES during
# the check process, as these PROAST variables are not used by MutSeqR

utils::globalVariables(
  names = c(
    # # PROAST internal functions
    # "f.adjust.saved", "f.cat", "f.ced.cat", "f.CI.sel.ord", "f.delete.gw",
    # "f.dtype6.mn", "f.expect.cat", "f.firstcheck", "f.nested.con",
    # "f.plot.cxt", "f.plot.sep", "f.press.key.to.continue", "f.quick.cat",
    # "f.select.m46.con", "f.show.cat", "f.start.cat", "f.use.comma.for.decimal",

    # # PROAST internal state variables
    # ".gw.size", ".Pr.last", ".ypos"
    "priority_order", "variation_type", "total_depth", ".N", ".", "contig", 
    "proportion_min", "proportion_max", "CED", "Selected.Model", "CEDL", "CEDU", 
    "weights", "dum", "vaf", "filter_reason", "in_regions", "x.circle_coords", 
    "y.circle_coords", "id", "color_column", "label", "y", "radius", "Significance",
    "ref_level", "Symbol", "subtype", "subtype_class", "subtype_labels", "proportion",
    "Group", "ProportionPlot", "Samples", "mutation", "xmin", "xmax", "ymin", "ymax",
    "x_variable"
  )
)