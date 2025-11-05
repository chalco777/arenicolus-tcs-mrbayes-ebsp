from Bio import Entrez, SeqIO

codes = ["SGR3","SGR4","MVZ241596","MVZ237413","JWA338","JWA470","CAS229140","CAS223822","MVZ149956"]


def search_accessions(gene_terms):
    results = {}
    for code in codes:
        query = (
            f'"Sceloporus graciosus"[Organism] AND ({gene_terms}) AND ({code}[All Fields]) '
            f'AND (KC853700:KC854199[ACCN])'
        )        
        handle = Entrez.esearch(db="nuccore", term=query, retmax=200)
        ids = Entrez.read(handle)["IdList"]; handle.close()
        for gi in ids:
            print (f"Fetching GI {gi} for code {code}...")
            gb = Entrez.efetch(db="nuccore", id=gi, rettype="gb", retmode="text")
            rec = SeqIO.read(gb, "genbank"); gb.close()
            acc = rec.annotations.get("accessions", [rec.id])[0]
            print(rec)

            # We get the description
            defline = rec.description.lower()
            print(defline)

            is_cytb = ("cytochrome b" in defline) or ("cytb" in defline) or ("cob" in defline)
            is_nd1  = ("nd1" in defline) or ("nadh dehydrogenase subunit 1" in defline) or ("nadh1" in defline) or ("nadh" in defline)

            #only if the defline indicates
            if is_cytb:
                gene = "cytb"
                results.setdefault(gene, {}).setdefault(code, set()).add(acc)
            elif is_nd1:
                gene = "nd1"
                results.setdefault(gene, {}).setdefault(code, set()).add(acc)
            # si no matchea, no añadimos nada
            else:
                print(f"Warning: could not determine gene for GI {gi} with defline: {defline}")
    return results

res = {}
res.update(search_accessions("cytochrome b[Title] OR cytb[Title]"))
res.update(search_accessions('ND1[Gene] OR "NADH dehydrogenase subunit 1"[Title] OR NADH1[Title]'))

for gene in ("cytb","nd1"):
    fn = f"../data/sequences_by_gene/acc_outgroups_{gene}.tsv"
    with open(fn, "w") as fh:
        for code, accs in sorted(res.get(gene, {}).items()):
            for acc in sorted(accs):
                fh.write(f"{code}\t{acc}\n")
    print("Written:", fn)
