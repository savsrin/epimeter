#!/usr/bin/env python3
"""
PeptideIndex
@author: savsr
PeptideIndex allows users to build a PeptideIndex object
that stores AnnoyIndexes. Each AnnoyIndex contains peptides
of a specific k-mer length processed from a protein FASTA 
file. 
"""
#change to order of blosum matrix 
import os
import errno
import sys
from annoy import AnnoyIndex
from collections import defaultdict  
class PeptideIndex(object):
    _amino_acids = {'A': 0, 'R': 1, 'N' : 2, 'D' : 3, 'C' : 4, 'Q' : 5, 
               'E' : 6, 'G' : 7, 'H' : 8, 'I' : 9, 'L' : 10, 'K' : 11,
               'M': 12, 'F' : 13, 'P': 14, 'S' : 15,'T' : 16, 'W' : 17,
               'Y': 18,'V' : 19, 'B' : 20, 'Z' : 21, 'X' : 22, '*': 23, 'U': 24} 
               #U is selenocystine; since it is very rare in human proteins and 
               #has a chemical structure similar to cysteine it will have the 
               #same scores assigned in the default Blosum matrix 

    _amino_acids_string = "ARNDCQEGHILKMFPSTWYVBZX*U"

    #default scores matrix with blosum 62 scores
    _scores =           [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4, 0],
                        [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4, -3],
                        [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4, -3], 
                        [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4, -3],  
                        [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4, 9],
                        [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4, -3], 
                        [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4, -4],
                        [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4, -3],
                        [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4,-3],
                        [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4, -1],
                        [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4, -1],
                        [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4, -3],
                        [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4,-1],
                        [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4,-2],
                        [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4, -3],
                        [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4, -1],
                        [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4, -1],
                        [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4, -2],
                        [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4, -2],
                        [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4, -1],
                        [-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4, -3],
                        [-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4, -3],
                        [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4, -2],
                        [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4, -4],
                        [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4, 9]]

    _weights = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0]#default value
    
    def __init__(self): 
        """"initializes PeptideIndex with dictionary that stores AnnoyIndexes,
            a directory name where the PeptideIndex will be saved, a boolean
            variable that keeps track of whether that index has been saved, and
            the default values for the scores matrix and weights vector 
        """
        self.annoy_indexes = {}
        self.index_dir = None
        self.n_trees = 100
        self.saved = False 
        self.peptide_nns = ()
        self.scores = self._scores
        self.weights = self._weights 
        
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
                                          self._get_vector(peptide)) 
        except KeyError:
            self.annoy_indexes[k] = AnnoyIndex(k,
                               metric='blosum', 
                               num_amino_acids=len(self._amino_acids_string), 
                               weights=self.weights,
                               scores=self.scores)
            self.annoy_indexes[k].add_item(item_number,
                                          self._get_vector(peptide))
    
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


    def set_scores_and_weights(self, scores_file, weights_file=""):
        """ sets the scores matric and weights vector to user
            specified values 
            note: parsing of weights_file and resetting self.weights
            has not been implemented yet 

            scores_file: path to file with user scores matrix
            weights_file: path to file with user weights vector

            Return Value: none 
        """
        if scores_file is not None: 
            #reinitialize the default amino acid string 
            #and default scores matrix 
            self._amino_acids_string = ""
            self.scores = [] 
            #temporary dictionary to store amino acids from user
            amino_acids = defaultdict(int) 
            #parse file with user scores matrix 
            with open(scores_file) as scores_file:
                read_amino_acids = False
                for line in scores_file:
                    #skip comments in file
                    if line.startswith("#"):
                       next(scores_file)
                    elif read_amino_acids == False: 
                        # header contains amino acids from matrix 
                        header = line.split() 
                        aa_count = 0
                        for char in header:
                            #set amino acid dict and string
                            #to have amino acids from user's matrix
                            self._amino_acids_string +=  char
                            amino_acids[char] = aa_count
                            aa_count = aa_count + 1 
                        read_amino_acids = True
                    else: 
                        header = (line.strip()[1:]).split()
                        num_amino_acids= len(header)
                        row = []
                        for i in range(num_amino_acids):
                            row.append(float(int(header[i])))
                        self.scores.append(row) 
            self._amino_acids = amino_acids 
         
        #print(self.scores)
          
    def save(self, index_name):
        """saves AnnoyIndexes to a directory named index_name and the files
           with the corresponding scores matrix and weight vector 
           (note: only scores has been implemented so far)
            
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
        for key in self.annoy_indexes:
            annoy_index_file = os.path.join(self.index_dir, str(key) + ".pid")  
            self.annoy_indexes[key].save(annoy_index_file)
        #writes and saves a file with the blosum matrix in the directory with
        #the annoy indices; blosum matrix format contains header with amino acids
        #and then before each row of scores is the corresponding amino acid
        scores_file = open(os.path.join(self.index_dir,"scores.txt"), "w")
        #writes header with amino acids
        scores_file.write("  ") 
        for i in range(len(self._amino_acids_string)):
            scores_file.write(self._amino_acids_string[i])
            if i < (len(self._amino_acids_string) - 1):
                scores_file.write(" ")
        scores_file.write("\n")
        #writes the rows of scores with the corresponding amino acid 
        for row in range(len(self.scores)):
            scores_file.write(self._amino_acids_string[row])
            for num in range(len(self.scores[row])):
                scores_file.write(" " + str(int(self.scores[row][num])))
            scores_file.write("\n")
        scores_file.close() 

        self.saved = True
                      
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
   
    def get_nns_by_epitope(self, epitope, num_neighbors=100, 
                           search_k=-1, include_distances=True):
        """ queries an epitope in appropriate AnnoyIndex and returns nearest 
            neighbors 
            
            epitope: epitope (string) being queried for neighbors
            num_neighbors: number of epitope matches found, default = 8
            search_k: number of nodes AnnoyIndex inspects while searching
            include_distances: determines whether distances between query
            and neighbors are returned along with the nearest neighbors
            
            Return Value: tuple, contains nearest neighbors, distances, and
            original peptides for the nearest neighbors
        """
        k = len(epitope)
        results = ()
        try: 
            temp_results = (self.annoy_indexes[k].get_nns_by_vector(
                                                 self._get_vector(epitope),
                                                 n=num_neighbors, 
                                                 search_k=search_k, 
                                                 include_distances
                                                 =include_distances))
        except KeyError: 
            a = self._lazy_load(str(k))
            #if blosum file exists, set the blosum matrix here 
            temp_results =  (a.get_nns_by_vector(self._get_vector(epitope),
                                          n=num_neighbors, 
                                          search_k=search_k, 
                                        include_distances=include_distances))

        results = results + (temp_results[0],) #adds item numbers for nns
        results = results + (temp_results[1],) #adds distances for nns
        
        peptides = [self._get_peptide(k, item_num) 
                   for item_num in results[0]]
        results = results + (peptides,) #adds peptide sequences for nns
        return results 


    def get_peptide_nns(self): 
        return self.peptide_nns
    
    def _get_peptide(self, k, item_num): 
        """provides corresponding peptide for item number in AnnoyIndex
           
           k: length of kmer
           item_num: item number of peptide vector in the annoy index
           
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
        scores_file_name = os.path.join(self.index_dir,"scores.txt")
        #resets scores/weights for the loaded annoy index 
        if os.path.exists(scores_file_name):
            self.set_scores_and_weights(scores_file_name, "")
        a = AnnoyIndex(int(k), metric='blosum',
                       num_amino_acids=len(self._amino_acids_string), 
                       weights=self.weights, 
                       scores=self.scores)
        #throws FileNotFound exception if file not present
        a.load(os.path.join(self.index_dir, k + ".pid"))
        self.annoy_indexes[int(k)] = a
        return a 
