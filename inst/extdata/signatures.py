try:
    from SigProfilerAssignment import Analyzer as Analyze
except Exception as e:
    print(f"Importing SigProfilerAssignment.Analyzer FAILED: {e}")
try:
    from SigProfilerExtractor import sigpro as sig
except Exception as e:
    print(f"Importing SigProfilerExtractor.sigpro FAILED: {e}")
try:
    from SigProfilerMatrixGenerator import install as genInstall
except Exception as e:
    print(f"Importing SigProfilerMatrixGenerator.install FAILED: {e}")
try:
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
except Exception as e:
    print(f"Importing SigProfilerMatrixGenerator.scripts.SigProfilerMatrixGeneratorFunc FAILED: {e}")

def install_MutSeqR(genome,
                    custom=False,
                    rsync=False,
                    bash=True,
                    ftp=True,
                    fastaPath=None,
                    transcriptPath=None,
                    exomePath=None,
                    offline_files_path=None,
                    volume=None):
    genInstall.install(genome = genome,
                       custom=custom,
                       rsync=rsync,
                       bash=bash,
                       ftp=ftp,
                       fastaPath=fastaPath,
                       transcriptPath=transcriptPath,
                       exomePath=exomePath,
                       offline_files_path=offline_files_path,
                       volume=volume)

def matrix_generator_MutSeqR(project,
                    reference_genome,
                    path_to_input_files,
                    plot = True,
                    exome = False,
                    bed_file = None,
                    chrom_based = False,
                    tsb_stat = True,
                    seqInfo = True,
                    cushion = 100):
   matGen.SigProfilerMatrixGeneratorFunc(project = project,
                    reference_genome = reference_genome,
                    path_to_input_files = path_to_input_files,
                    plot = plot,
                    exome = exome,
                    bed_file = bed_file,
                    chrom_based = chrom_based,
                    tsb_stat = tsb_stat,
                    seqInfo = seqInfo,
                    cushion = cushion)

def cosmic_fit_MutSeqR(samples,
                       output,
                       input_type="matrix",
                       context_type="96",
                       cosmic_version=3.4,
                       exome=False,
                       genome_build="GRCh38",
                       signature_database=None,
                       exclude_signature_subgroups=None,
                       export_probabilities=True,
                       export_probabilities_per_mutation=False,
                       make_plots=True,
                       sample_reconstruction_plots = "png",
                       verbose=True):
  Analyze.cosmic_fit(samples,
                       output,
                       input_type=input_type,
                       context_type=context_type,
                       cosmic_version=cosmic_version,
                       exome = exome,
                       genome_build=genome_build,
                       signature_database=signature_database,
                       exclude_signature_subgroups=exclude_signature_subgroups,
                       export_probabilities=export_probabilities,
                       export_probabilities_per_mutation=export_probabilities_per_mutation,
                       make_plots=make_plots,
                       sample_reconstruction_plots = sample_reconstruction_plots,
                       verbose=verbose
                       )
print("Defined cosmic_fit_MutSeqR")

def denovo_fit_MutSeqR(samples,
                       output,
                       input_type="matrix",
                       signatures="path/to/input/denovo/signatures/file",
                       genome_build="GRCh38"):
    Analyze.denovo_fit(samples,
                        output,
                        input_type=input_type,
                        signatures=signatures,
                        genome_build=genome_build)
print("Defined denovo_fit_MutSeqR")

def decompose_fit_MutSeqR(samples,
                       output,
                       input_type="matrix",
                       signatures="path/to/input/denovo/signatures/file",
                       signature_database="path/to/optional/reference/signatures/database/file",
                       genome_build="GRCh38"):
    Analyze.decompose_fit(samples,
                        output,
                        input_type=input_type,
                        signatures=signatures,
                        signature_database=signature_database,
                        genome_build=genome_build)
print("Defined decompose_fit_MutSeqR")
