#' @title Set up Python environment for MutSeqR
#'
#' @description This function initializes the Python environment used by
#' MutSeqR.
#' It is not run automatically to avoid issues during installation and
#' checks.
#' @param force Logical. Whether to force reconfiguration even if an
#' environment already exists.
#' @returns None
setup_mutseqr_python <- function(force = FALSE) {
  reticulate::configure_environment("MutSeqR", force = force)
}

.onAttach <- function(libname, pkgname) {
  print_ascii_art()
}
#' This function prints ASCII art when the package is loaded
#' @returns None
print_ascii_art <- function() {
    art_path <- system.file("extdata", "ASCII_art_MutSeqR.txt", package = "MutSeqR")
    art <- readLines(art_path)
    packageStartupMessage(paste(art, collapse = "\n"))
}

#' Get Mutation Palette
#'
#' Internal helper function to determine the color palette based on resolution.
#' Supports types, base_6, base_12, base_96, base_192
#' @keywords internal
get_mutation_palette <- function(custom_palette = NULL,
                                 subtype_resolution,
                                 num_colours) {
    if (!is.null(custom_palette)) {
        # Check if custom_palette is a string matching a Brewer palette name
        is_brewer_name <- is.character(custom_palette) &&
            length(custom_palette) == 1 &&
            (custom_palette %in% rownames(RColorBrewer::brewer.pal.info))

        if (is_brewer_name) {
            max_colors <- RColorBrewer::brewer.pal.info[custom_palette, "maxcolors"]
            if (num_colours > max_colors) {
                grDevices::colorRampPalette(
                    RColorBrewer::brewer.pal(n = max_colors, name = custom_palette)
                )(num_colours)
            } else {
                p <- RColorBrewer::brewer.pal(n = num_colours, name = custom_palette)
                return(p)
            }
        } else {
            return(custom_palette)
        }
    } else {
        palette_types <- RColorBrewer::brewer.pal(8, "BrBG")
        names(palette_types) <- c(
            "complex", "deletion", "insertion", "mnv", "snv",
            "sv", "ambiguous", "uncategorized"
        )

        # Base 6 Roots (Spectral)
        # Order: T>G, T>C, T>A, C>T, C>G, C>A
        root_6 <- RColorBrewer::brewer.pal(6, "Spectral")
        names(root_6) <- c("T>G", "T>C", "T>A", "C>T", "C>G", "C>A")

        # Base 12 Roots (Sanger: Reds, Greys, Blues, Greens)
        # Order: T>G, T>C, T>A (Reds) | G>T, G>C, G>A (Greys) | ...
        root_12 <- c(
            RColorBrewer::brewer.pal(3, "Reds"),
            RColorBrewer::brewer.pal(3, "Greys"),
            RColorBrewer::brewer.pal(3, "Blues"),
            RColorBrewer::brewer.pal(3, "Greens")
        )
        names(root_12) <- c(
            "T>G", "T>C", "T>A", "G>T", "G>C", "G>A",
            "C>T", "C>G", "C>A", "A>T", "A>G", "A>C"
        )

        # This function takes a vector of root colors and expands each one into 'n'
        # shades. It creates a gradient from a lighter version to a darker version
        # of the root.
        expand_colors <- function(roots, n = 16, dict) {
            expanded <- unlist(lapply(roots, function(hex_color) {
                light_end <- grDevices::colorRampPalette(c("white", hex_color))(10)[4]
                dark_end <- grDevices::colorRampPalette(c(hex_color, "black"))(10)[4]
                grDevices::colorRampPalette(c(light_end, dark_end))(n)
            }))
            names(expanded) <- dict
            return(expanded)
        }

        palette_lookup <- list(
            "none" = "black",
            "type" = palette_types,
            "base_6" = c(root_6, palette_types),
            "base_12" = c(root_12, palette_types),
            "base_96" = c(expand_colors(root_6, n = 16, dict = MutSeqR::subtype_list[[subtype_resolution]]), palette_types),
            "base_192" = c(expand_colors(root_12, n = 16, dict = MutSeqR::subtype_list[[subtype_resolution]]), palette_types)
        )
        return(palette_lookup[[subtype_resolution]])
    }
}
