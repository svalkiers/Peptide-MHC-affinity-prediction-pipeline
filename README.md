# Peptide-MHC affinity prediction pipeline (PMAPP)

To be able to run the PMAPP, you need to install [netMHCpan4.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan) and [netCTLpan](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netCTLpan). Make sure that you add the models to your `$PATH` variable.

Download the peptide-MHC_affinity_prediction.py script from the code folder. Open the script in your Python interpreter, edit the parameters of the function `pipeline(model, path, allele_file, FASTA_file)` and run the script.

Arguments:
- `model`: The model you want to use. Right now, the pipeline supports *'netMHCpan'* and *'netCTLpan'*.
- `path`: The path to your current working directory.
- `allele_file`: A file with HLA alleles. You can type `netMHCpan -listMHC` or `netCTLpan -listMHC` in your terminal to check the available HLA alleles.
- `FASTA_file`: A FASTA file of the amino acid sequences you want to screen.

Alternatively, you can run the script in your terminal by using the `python3 peptide-MHC_affinity_prediction.py` command.

*Note*: This script will only work on a **Linux system**.
