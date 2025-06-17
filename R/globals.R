# This file tells R CMD check about variables used in the context
# of non-standard evaluation
utils::globalVariables(
  names = c(
    "priority_order", "variation_type", "total_depth", ".N", ".", "contig", 
    "proportion_min", "proportion_max", "CED", "Selected.Model", "CEDL", "CEDU", 
    "weights", "dum", "vaf", "filter_reason", "in_regions", "x.circle_coords", 
    "y.circle_coords", "id", "color_column", "label", "y", "radius", "Significance",
    "ref_level", "Symbol", "subtype", "subtype_class", "subtype_labels", "proportion",
    "Group", "ProportionPlot", "Samples", "mutation", "xmin", "xmax", "ymin", "ymax",
    "x_variable"
  )
)