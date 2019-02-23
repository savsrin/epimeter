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
from annoy import AnnoyIndex
class PeptideIndex(object):
    _amino_acids = {'A': 0, 'B': 1, 'C' : 2, 'D' : 3, 'E' : 4, 'F' : 5, 
               'G' : 6, 'H' : 7, 'I' : 8, 'K' : 9, 'L' : 10, 'M' : 11,
               'N': 12, 'P' : 13, 'Q': 14, 'R' : 15,'S' : 16, 'T' : 17,
               'U': 18,'V' : 19, 'W' : 20, 'X' : 21, 'Y' : 22, 'Z': 23}
    _amino_acids_string = "ABCDEFGHIKLMNPQRSTUVWXYZ"
    
    def __init__(self): 
        """"initializes PeptideIndex with dictionary that stores AnnoyIndexes,
            a directory name where the PeptideIndex will be saved, and boolean
            variable that keeps track of whether that index has been built
        """
        self.annoy_indexes = {}
        self.index_dir = ""
        self.n_trees = 10
        self.saved = False 
        
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
            self.annoy_indexes[k].add_item(item_number, 
                                          self.__get_vector(peptide)) 
        except KeyError:
            self.annoy_indexes[k] = AnnoyIndex(k,
                               metric = 'euclidean')
            self.annoy_indexes[k].add_item(item_number,
                                          self.__get_vector(peptide))
    
    def build(self):
        """builds Peptide Index""" 
        for key in self.annoy_indexes:
            self.annoy_indexes[key].build(self.n_trees)
    
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
        if not os.path.exists(self.index_dir):
            os.makedirs(self.index_dir)
        self.build() # only built indexes are saved  
        for key in self.annoy_indexes:
            annoy_index_file = os.path.join(self.index_dir, str(key) + ".pid")  
            self.annoy_indexes[key].save(annoy_index_file)
        self.saved = True
                      
    def load(self, index_name): 
        """loads PeptideIndex if not built reinitializes annoy_indexes if 
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
        if self.index_dir != "": 
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
            
            Return Value: tuple, contains nearest neighbors and distances
        """
        k = len(epitope)
        try: 
            results = self.annoy_indexes[k].get_nns_by_vector(self.__get_vector(epitope),
                                                 n=num_neighbors, 
                                                 search_k=search_k, 
                                                 include_distances
                                                 =include_distances) 
        except KeyError: 
            a = self.__lazy_load(str(k))
            results = a.get_nns_by_vector(self.__get_vector(epitope),
                                          n=num_neighbors, 
                                          search_k=search_k, 
                                          include_distances=include_distances)
        return results 
    
    def __get_peptide(self, k, item_num): 
        """provides corresponding peptide for item number in AnnoyIndex;
           
           k = length of kmer
           i = item number of peptide vector in the annoy index
           
           Return Value: string
        """
        vector = self.annoy_indexes[k].get_item_vector(item_num)
        peptide = ""
        for num in vector: 
            peptide = peptide + self._amino_acids_string[num]
        return peptide
    
    def __get_vector(self, peptide): 
        """convert amino acids in peptide to numerical vector"""
        return [self._amino_acids[peptide[i]] for i in range(len(peptide))]    
    
    def __lazy_load(self, k):
        """loads annoy indexes as needed
           
            k: kmer length 
            
            Return Value: AnnoyIndex"""
        a = AnnoyIndex(int(k), metric='euclidean') 
        #throws FileNotFound exception if file not present
        a.load(os.path.join(self.index_dir, k + ".pid"))
        return a 
        
    

        
        
        
    
        
        
            
            