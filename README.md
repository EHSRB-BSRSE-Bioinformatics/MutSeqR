# Duplex Sequencing Analysis

## Pre-abmle

To do: write a bit about the project, the scope, why it's necesary.

Provide a little background on duplex sequencing, describe the technology
briefly, and give some context for the type of data we are meant to be
processing with this package.

## Data import

The first step is importing data. We have written a function called 
`import_mut_data` to accomplish this in one line. Since the main goal of this 
package is to generate summary statistics, visualization, exploratory analysis, 
and other post-processing tasks such as mutational signature analysis or 
generalized linear modeling, the main piece of information you want to import is
the `.mut` file, the schema for which is described here:

(example mut file with columns and their descriptions?)

It is critical that this .mut file gets annotated with information such as the
genomic region from which the mutation originates, some information about the 
sample from which the mutation orignates, and that each row (mutation) recieves
a calculated group depth and group frequency (of interest for later analyses).

We do this by first importing the `.mut` file as a data frame, and then joining 
it with that other data (e.g., regions file), and finally converting this to a 
`granges` object. This facilitates use in other packages and makes doing 'genome
 math' on the ranges significantly easier.

### Metadata: an important consideration

The other important component of importing your data for proper use is to assign
 each mutation to a biological sample, and also make sure that some additional 
 information about each sample is present (e.g., a chemical treatment, a dose, 
 etc.). This is done by providing a sample data file (tab delimited, comma 
 delimited, etc.; the choice is up to the user, but the delimiter of the file 
 must be specified as a parameter in the function). Importantly, this is a file
 that would be analogous to "colData", or "column data", a term often used in 
 the `DESeq2` package. Hence, it must contain some information about an existing
 column in your `.mut` file, which is typically going to be sample. So the first
 column in your sample data file should indeed be `sample`. Then, additional 
 columns such as `dose` or `tissue` or `treatment` can be added, and these 
 columns will be joined with your `.mut` file to capture that information and
  associate it with each mutation.
  
### Other Notes

When first importing a `.mut` file, it is preferred to keep non-variant rows.
This allows the calculuation of mutation frequencies. The data set can be pared
down later to include only mutations of interest (SNVs, indels, SVs, or any 
combination).

To be filled in more...

## Step 2

## Step 3