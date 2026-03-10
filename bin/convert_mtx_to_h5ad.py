import argparse
import scanpy as sc

def parse_args():

    parser = argparse.ArgumentParser(description="Parse your predictions for evaluation.")

    # Single path specifying were mtx.gz data is deposited
    parser.add_argument(
        "--countsDirectory",
        required=True,
        type=str,
        help="Path to mtx data in gz format"
    )

    # Sample name
    parser.add_argument(
        "--sampleID",
        type=str,
        default="Sample",
        help="Sample name for output file"
    )

    return parser.parse_args()


def main():
    # Get arguments
    args = parse_args()

    # Read data
    print(f"Directory: {args.countsDirectory}")
    adata = sc.read_10x_mtx(args.countsDirectory, var_names="gene_ids", cache=True)

    # Annotate transcripts other than CDSs
    # adata.var["tRNA"] = adata.var['gene_ids'].str.contains("_tRNA_")
    # adata.var["rRNA"] = adata.var['gene_ids'].str.contains("_rRNA_")
    # adata.var["other_ncRNA"] = adata.var['gene_ids'].str.contains("_RNA_")

    # Write Anndata to output
    adata.write_h5ad(f"{args.sampleID}_anndata.h5ad")


if __name__ == '__main__':
    main()
    