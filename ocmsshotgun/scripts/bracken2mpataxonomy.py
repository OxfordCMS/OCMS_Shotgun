
from argparse import ArgumentParser
import subprocess
import os
import re

def main(argv=None):
    '''
    translates NCBI taxids in default kraken2 output to names in 'mpa' style using Taxonkit
    this script is modified from CGATMetaSequencing written by Matt Jackson
    https://github.com/microbialman/CGATMetaSequencing/blob/master/scripts/translateKraken2.py
    '''
    
    #get the default kraken2 output and the file to save translated names file to
    
    parser = ArgumentParser()
    parser.add_argument("--mergedbracken",dest="infile", help="merged bracken output file to be translated")
    parser.add_argument("--translatedout", dest="outfile", help="Output file name")
    parser.add_argument("--taxdatadir", dest="tddir", help="Directory containing names.dmp and nodes.dmp from NCBI taoxnomy (for TaxonKit, uses taxonkit default if not set).", )
    parser.add_argument("--taxaprefixes", dest="tpre", help="Set to True is the taxonomy databases already includes p__ etc. in its taxonomy names (e.g. GTDB databases).")
    
    args = parser.parse_args()

    #open the files
    infile = open(args.infile, 'r')
    outfile = open(args.outfile, 'w')
        
    #read name to taxid mapping
    readnames=[]
    taxids=[]

    #unique taxids to find the taxonomic names for
    uniqueids=[]

    #go through merged bracken result and get taxids
    pattern = r"(.*)-(\d+)$"

    for n,i in enumerate(infile):
        row = i.split("\t")
        if row[0] == 'bin':
            continue
        match = re.search(pattern, row[0])
        if match:
            taxid = match.group(2)
            readname = "%s taxid %s" % (match.group(1), taxid)

            #readname = row[1]
            #taxid = row[2]
            readnames.append(readname)
            taxids.append(taxid)
            if taxid not in uniqueids:
                uniqueids.append(taxid)
    infile.close()
    
    #dictionary to store found ids
    taxdic={}

    #write the taxids to a temp file
    tempfilename=os.getcwd()+"/"+args.outfile+"_TEMP_Taxids.txt"
    tempfile=open(tempfilename,"w")
    tempfile.write("\n".join(uniqueids))
    tempfile.close()

    #get lineage names using TaxonKit
    if args.tddir:
        command="cat {} | taxonkit lineage --data-dir {} | taxonkit reformat --data-dir {}".format(tempfilename, args.tddir, args.tddir)
    else:
        command="cat {} | taxonkit lineage | taxonkit reformat".format(tempfilename)
    p=subprocess.check_output(command, shell=True)
    taxonkit=p.decode().split("\n")

    #delete tempfile
    subprocess.run(["rm",tempfilename])

    #taxonomic levels
    levs=["k__","p__","c__","o__","f__","g__","s__"]

    #function to generate the mpa name from a given taxid by calling Taxonkit
    def formName(name):
        names = name.split(";")
        formnames = []
        for i in range(len(names)):
            if names[i] != "":
                namestr=names[i]
                if args.tpre == "True":
                    namestr=namestr[3:]
                formnames.append(levs[i]+namestr.replace(" ","_"))
            else:
                formnames.append(levs[i]+"unassigned")
        return(";".join(formnames))

    #reformat to mpa style
    for i in taxonkit:
        if i != "":
            row=i.split("\t")
            tid=row[0]
            mpaname=formName(row[2])
            if mpaname != "":
                taxdic[tid]=mpaname

    #write the translated read to taxname file
    for i in range(len(readnames)):
        t=taxids[i]
        if t in taxdic:
            outfile.write("{}\t{}\n".format(readnames[i],taxdic[t]))

    if argv is None:
        argv = sys.argv

if __name__ == "__main__":
    sys.exit(main(sys.argv))
