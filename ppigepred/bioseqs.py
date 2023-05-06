import pyensembl


def get_ensembl_release(assembly):
    ensembl = pyensembl.reference_name.genome_for_reference_name(assembly)
    ensembl.download()
    ensembl.index()
    return ensembl


def gene_ids_to_protein_ids(ensembl, gene_ids):
    genes = (ensembl.gene_by_id(id) for id in gene_ids)
    return {transcript.protein_id
            for gene in genes
            for transcript in gene.transcripts
            if transcript.protein_id}


def protein_ids_to_gene_ids(ensembl, protein_ids):
    return {prot_id: ensembl.gene_id_of_protein_id(prot_id)
            for prot_id in protein_ids}