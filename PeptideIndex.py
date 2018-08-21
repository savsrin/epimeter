#!/usr/bin/env python3
"""
PeptideIndex
@author: savsr
PeptideIndex allows users to build a PeptideIndex object
that stores AnnoyIndexes. Each AnnoyIndex contains peptides
of a specific k-mer length processed from a protein FASTA 
file. 
"""

import os
import errno
from annoy import AnnoyIndex
class PeptideIndex(object):
    _amino_acids = {'A': 0, 'R': 1, 'N' : 2, 'D' : 3, 'C' : 4, 'Q' : 5, 
               'E' : 6, 'G' : 7, 'H' : 8, 'I' : 9, 'L' : 10, 'K' : 11,
               'M': 12, 'F' : 13, 'P': 14, 'S' : 15,'T' : 16, 'W' : 17,
               'Y': 18,'V' : 19, 'B' : 20, 'Z' : 21, 'X' : 22, '*': 23, 'U': 24}
    _amino_acids_string = "ARNDCQEGHILKMFPSTWYVBZX*"
    
    def __init__(self): 
        """"initializes PeptideIndex with dictionary that stores AnnoyIndexes,
            a directory name where the PeptideIndex will be saved, and boolean
            variable that keeps track of whether that index has been saved
        """
        self.annoy_indexes = {}
        self.index_dir = None
        self.n_trees = 150
        self.saved = False 
        #self.peptide_neighbors = []
        
    def add_item(self, item_number, peptide):
        """adds peptides to the appropriate AnnoyIndex in the PeptideIndex
        
           item_number:specificies which peptide # is being added to the index 
           peptide: protein kmer that is being indexed
           
           Return value: none
        """
        if self.saved:
            # annoy raises an error that vector can't be modified
            raise ValueError("PeptideIndex already saved, can't add item.")
        k = len(peptide)   #kmer length
        try:
            v = self._get_vector(peptide)
            self.annoy_indexes[k].add_item(item_number, 
                                          v) 
            del v 
        except KeyError:
            self.annoy_indexes[k] = AnnoyIndex(k,
                               metric='blosum')
            v = self._get_vector(peptide)
            self.annoy_indexes[k].add_item(item_number,
                                          v)
            del v 

    
    def build(self):
        """builds Peptide Index""" 
        print("starting to build (in peptide class) \n")
        for key in self.annoy_indexes:
            self.annoy_indexes[key].build(self.n_trees)
        print("finished building in peptide class \n")
    
    def set_n_trees(self, n_trees):
        """ sets trees to user specified value"""
        if self.saved: 
            raise ValueError("PeptideIndex has already been saved,"
                             + "can't set the number of trees.")
        self.n_trees = n_trees
                        
    def save(self, index_name):
        """saves AnnoyIndexes to a directory named index_name
            
        index_name: path to directory
        
        Return Value: none
        """
        self.index_dir = index_name
        try:
            os.makedirs(self.index_dir)
        except OSError as e:
            if e.errno == errno.EEXIST and os.path.isdir(self.index_dir):
                pass
            else:
                raise
        self.build() # only built indexes are saved  
        print("starting to save indices \n")
        for key in self.annoy_indexes:
            annoy_index_file = os.path.join(self.index_dir, str(key) + ".pid")  
            self.annoy_indexes[key].save(annoy_index_file)
        self.saved = True
        print("finished saving indices (in peptide class)")
                      
    def load(self, index_name="./"): 
        """loads PeptideIndex ; reinitializes annoy_indexes if 
           directory name has been changed
        
           index_name: path to directory where PeptideIndex is stored
           
           Return Value: none 
        """
        #raise an error, because loading into this PeptideIndex
        # will destroy this object.
        #items may have been added to this index so it may have not
        #been built and saved. 
        if not self.saved and any(self.annoy_indexes):  
            raise ValueError("PeptideIndex has not been saved.") 
        if self.index_dir is not None: 
            #reinitialize the PeptideIndex because we may be loading
            # a different peptide index into this object (eg., if load
            # is called twice)
            self.annoy_indexes = {} 
        self.index_dir = index_name 
        self.saved = True 
   
    def get_nns_by_epitope(self, epitope, num_neighbors=8, 
                           search_k=-1, include_distances=True):
        """ queries an epitope in appropriate AnnoyIndex and returns nearest 
            neighbors 
            
            epitope: epitope (string) being queried for neighbors
            num_neighbors: number of epitope matches found, default = 8
            search_k: number of nodes AnnoyIndex inspects while searching
            include_distances: determines whether distances between query
            and neighbors are returned along with the nearest neighbors
            
            Return Value: tuple, contains nearest neighbors, distances, and 
            peptide sequences for nearest neighbors 
        """
        k = len(epitope)
        results = () 
        try: 
            temp_results = self.annoy_indexes[k].get_nns_by_vector(
                                                 self._get_vector(epitope),
                                                 n=num_neighbors, 
                                                 search_k=search_k, 
                                                 include_distances
                                                 =include_distances) 
        except KeyError: 
            a = self._lazy_load(str(k))
            temp_results = a.get_nns_by_vector(self._get_vector(epitope),
                                          n=num_neighbors, 
                                          search_k=search_k, 
                                          include_distances=include_distances)
        
    
        results = results + (temp_results[0],) #adds item numbers for nns
        results = results + (temp_results[1],) #adds distances for nns 
        peptides = [self._get_peptide(k, item_num) 
                   for item_num in results[0]]
        #print(peptides)
        results = results + (peptides,) #adds peptide sequences for nns

        return results 
    
    def _get_peptide(self, k, item_num): 
        """provides corresponding peptide for item number in AnnoyIndex
           
           k: length of kmer
           i: item number of peptide vector in the annoy index
           
           Return Value: string
        """
        vector = self.annoy_indexes[k].get_item_vector(item_num)
        peptide = ""
        for num in vector: 
            peptide = peptide + self._amino_acids_string[int(num)]
        return peptide
    
    def _get_vector(self, peptide): 
        """convert amino acids in peptide to numerical vector"""
        return [self._amino_acids[peptide[i]] for i in range(len(peptide))]    
    
    def _lazy_load(self, k):
        """loads annoy indexes as needed
           
            k: kmer length as a string
            
            Return Value: AnnoyIndex
        """
        a = AnnoyIndex(int(k), metric='blosum') 
        #throws FileNotFound exception if file not present
        a.load(os.path.join(self.index_dir, k + ".pid"))
        self.annoy_indexes[int(k)] = a
        return a 
