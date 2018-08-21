#!/usr/bin/env python3
"""
epimeter
@author: savsr

epimeter is a tool that determines how similar a given neo-epitope is to the 
human proteome by calculating the distance between the peptide sequences. 
epimeter has two functions: "index" and "query." The "index" function takes 
a FASTA file of proteins as an input and indexes all of the protein k-mers
(user specifies k) into an Annoy Index. The "query" function then allows the
user to search for a specific neo-epitope candidate within the Annoy Index
allowing for a much more stream-lined and efficient search process than 
current bioinformatics search tools such as BLAST. 
"""
import PeptideIndex 
import argparse
from collections import defaultdict  
import sys #TODO remove

def build_index(protein_fasta, kmer_bounds, index_name):  
    """builds protein from protein fasta, indexes kmers
    in protein into a Peptide Index, and saves the Peptide Index.
    Each set of kmers has its own AnnoyIndex in the Peptide Index
    
    protein_fasta: path to file of protein seqs
    kmer_bounds: range of kmers for which AnnoyIndexes need to be built
    index_name: name of dir where PeptideIndex will be saved 
    
    Return Value: PeptideIndex 
    """
    print(kmer_bounds)
    min_k = int(kmer_bounds[0])
    max_k = int(kmer_bounds[1])
    peptide_index = PeptideIndex.PeptideIndex()
    protein = [] 
    with open(protein_fasta) as protein_fasta:
        #skip the first header in file
        next(protein_fasta)
        #keeps peptide count for each kmer AnnoyIndex in PeptideIndex
        item_numbers = defaultdict(int) 
        #for testing purposes: count = 0 
        for line in protein_fasta:   
            if line.startswith('>'):
                # New protein! Process k-mers from last protein assembled
                # Start with the full protein
                protein = ''.join(protein)
                for k in range(min_k, max_k + 1):
                    for i in range(len(protein) - k + 1):
                        item_numbers[k] = item_numbers[k] + 1
                        peptide_index.add_item(item_numbers[k], 
                                               protein[i:i+k]) 
                        #count = count + 1
                        #if count % 5000 == 0: 
                           # print(str (count) + "\n")   
                            #sys.stdout.flush()       
                protein = []
             #continue assembling the protein
            else: 
                protein.append(line.strip())
        #update the index for the last protein in the file 
        if len(protein) > 0: 
            protein = ''.join(protein)
            for k in range(min_k, max_k + 1):
                for i in range(len(protein) - k + 1):
                    item_numbers[k] =item_numbers[k] + 1
                    peptide_index.add_item(item_numbers[k], 
                                           protein[i:i+k])
   # print("finished adding all items to peptide index, count = " 
       #   + str(count) + "\n")
    peptide_index.save(index_name)
    #print("finished saving peptide index \n")
    return peptide_index 

def query_epitope(epitopes, index_name):
    """ takes each epitope from a file of epitopes, queries for the
    nearest neighbors, and then writes the nearest neighbors + distances
    for each epitope per line in a .csv file
    
    epitopes: path to file of epitopes; assumes one epitope per line
    index_name: path to directory where PeptideIndex will be loaded from
    
    Return Value: none
    """

    nearest_neighbors = open("nearest_neighbors.csv", "w")
    peptide_neighbors = open("peptide_nns.csv", "w")
    peptide_index = PeptideIndex.PeptideIndex()
    peptide_index.load(index_name) 
    with open(epitopes) as epitopes:
        for line in epitopes: 
            epitope = line.strip()
            sys.stdout.flush() 
            #returns a tuple with the neighbors and the distances
            results = peptide_index.get_nns_by_epitope(epitope,  
                                                      num_neighbors = 50, 
                                                      search_k=-1, 
                                                      include_distances=True)
            #stores neighbors and distances from tuple as lists
            neighbors = results[0]
            distances = results[1]
            peptides = results[2]
            for i in range(len(neighbors)):
                #insert commas between the neighbors and respective distances 
                #for the epitope
                nearest_neighbors.write(str(neighbors[i]) + "," + str(distances[i]))   
                peptide_neighbors.write(str(neighbors[i]) + "," + peptides[i])
                if i < len(neighbors) - 1: 
                   nearest_neighbors.write(",")
                   peptide_neighbors.write("\n")
            nearest_neighbors.write("\n")
            peptide_neighbors.write("\n")
    nearest_neighbors.close()
    peptide_neighbors.close() 
    return peptide_index 

def main(): 
    #creates the top level parser 
    parser = argparse.ArgumentParser(
            description = "Epimeter can index an input protein file "
                          "or query an epitope"
        ) 
    subparsers = parser.add_subparsers(dest='subparser_name')
    #create the parser for the "index" subcommand
    index_parser = subparsers.add_parser('index')
    index_parser.add_argument('-p', '--protein-fasta',
                              type=str, required=True,
                              help="path to protein fasta")
                             
    index_parser.add_argument('-k', '--kmer-bounds', nargs=2,
                              help="bounds on kmer sizes to index, e.g., "
                                   "\"8 11\"", 
                              default=["8", "11"])
                              
    index_parser.add_argument('-i', '--index-name', required=True, 
                              help="path to output index")
    #create the parser for the "query" subcommand
    query_parser = subparsers.add_parser('query')
    query_parser.add_argument('-i', '--index-name', required=True, 
                              help="path to output index")
    query_parser.add_argument('-e', '--epitopes', type=str,
                              required=True, 
                              help="path to list of epitopes")
    #arguments user types into command line
    args = parser.parse_args()
   
    if args.subparser_name == 'index':
           build_index(args.protein_fasta, 
                       args.kmer_bounds,
                       args.index_name) 
    elif args.subparser_name == 'query': 
          query_epitope(args.epitopes,
                        args.index_name)
  
if __name__ == '__main__':
    main()
    
      
 
    
