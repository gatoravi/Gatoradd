/* CROWN MRCA Version 1
 * Written by Avinash using Andre Wehe's tree library Feb 11 2011
 * program to add taxa of a specified family to a phylogenetic tree using the CROWN option.
 * usage - ./executable initial_tree_file leaves_families_file no_of_replicates opfile_name 
     initial tree file = file containing single newick string
     leaves_families_file = file with leaves to be added. Format = leaf1 anc1 anc2
     no_of_replicates = number of output replicates desired
     opfile_name = name of output file
 * note: crown  => cannot add leaves which have the same two ancestors  i.e leaf1 anc1 anc1 is not allowed
 */
 
#include "common.h"
#include "argument.h"
#include "tree.h"
#include "tree_IO.h"
#include "tree_traversal.h"
#include "tree_subtree_info.h"
#include "tree_LCA.h"
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/foreach.hpp>
#include <math.h>

#define MAX_ADD_LEAF 1000 /* Max Number of new leaves that can be added */ 
#define ADDED 1
#define NOTADDED 0
#define MAX_TREE_TAXA 2000 /* Max number of taxa in the input tree */

using namespace std;
int binarysearch(double *larray, int size, double key);

int main(int argc, char* argv[]) { 
  
  if(argc!=5){
    cout<<"Incorrect Arguments ! Exiting! \n\t usage - ./exec tree_file leaves_file replicates opfile\n";
    exit(1);
  }
 
  /* set seed of Rand number generator */
  srand(time(NULL));
  
  /* Read in command line arguments */
  char* treefile = argv[1];
  char* leaves_file = argv[2];  
  int replicates = atoi(argv[3]);
  char* opfile = argv[4];
  
  /* Read initial starting tree from treefile */
  std::ifstream ifs;
  ifs.open (treefile);
  aw::Tree initial_t;
  aw::idx2name initial_t_name;
  aw::idx2weight_double initial_t_weight;   
  if(ifs.good()){       
    if (!aw::stream2tree(ifs, initial_t, initial_t_name, initial_t_weight)){
      cout<<"Unable to read tree from file! Exiting!";
      exit(1);
    }
  } else {
    cout<<"Unable to open tree file! Exiting!\n";
    exit(1);
  }
  
  /* Open output file */
  ofstream ofs(opfile);
  if(!ofs.good()) {
    cout<<"Unable to open output file! Exiting!\n";
    exit(1);
  }
  
  /* Create a set of all families */
  set<string> family_set;
  
  /* Max number of nodes in final tree */
  int max_total_nodes = 2*(MAX_TREE_TAXA + MAX_ADD_LEAF -1);
  
  /* Arrays to store parents of nodes and translate node indices */  
  int*    initial_parent_array = new int[max_total_nodes];
  unsigned int* translate_index = new unsigned int[max_total_nodes];
  
  /* Initial number of leaves and internal nodes */
  int leafn = 0;
  int initial_currnodecount = 0;
    
  /* Map name of leaves to id */
  map<std::string, int> name2id_map;
  
  /* Traverse tree to create parents array */
  TREE_POSTORDER2(k,initial_t){
    unsigned int currnode = k.idx;
    if(initial_t.is_leaf(currnode)) {
      name2id_map[initial_t_name[currnode]] = currnode;
      leafn++;
    }    
    initial_currnodecount++;
    initial_parent_array[currnode] = k.parent;     
  } 
  cout<<"\nThe number of leaves in the initial tree is "<<leafn;
    
  /* Create an lca object to find lcas of nodes. */
  aw::LCA lca;
  lca.create(initial_t);
  
  /* Create an array to store lca values of individual nodes */
  int leaf_lca_array[MAX_ADD_LEAF];
  string leaves_array[MAX_ADD_LEAF];
  
 
  /*number of leaves to be added*/
  int leafcount = 0;  
  
  /* Read in leaves to be added from file */
  std::ifstream ifs2(leaves_file);
  string current_leaf;
  if(ifs2.is_open()){
    getline(ifs2, current_leaf);
    string taxa, anc1, anc2;
    pair<set<string>::iterator, bool> family_iterator;
    
    while(ifs2.good()) {      
    
      /* Parse leaf to be added and its two ancestors */
      istringstream iss(current_leaf);
      getline(iss, taxa, '\t');
      getline(iss, anc1, '\t');
      getline(iss, anc2, '\t');
      
      /* Add family of current leaf to family set */
      string curr_leaf_family = anc1 + " " + anc2;
      family_iterator = family_set.insert(curr_leaf_family);
           
      /* Get IDs of current leaf's ancestors*/
      int anc1_id = name2id_map[anc1];
      int anc2_id = name2id_map[anc2];
      
      /*Identify families with single taxa */
      if(anc1_id == anc2_id) {
        cout<<"\nLEAF "<<leafcount+1<<" has same 2 Ancestors. Incorrect Input Specification. Exiting\n";
        exit(1);
      } 
      
      /* Find MRCA of ancestors of curr leaf */
      int lca_id = lca.lca(anc1_id, anc2_id);
      //cout<<endl<<"MRCA of curr leaf anc is: "<<lca_id;
      
      /* Store leaf label */
      leaves_array[leafcount] = taxa;
      
      /* Store lca in array */
      leaf_lca_array[leafcount++] = lca_id;
      cout<<endl<<leafcount<<" anc1_id is: "<<anc1_id<<" anc2_id is: "<<anc2_id<<" lca "<<lca_id;
      
      /* Check if LCA is root */
      if(lca_id == 0) {
        cout<<"\nLeaf "<<leafcount+1<<"ERROR - LCA is root Exiting !\n";
        exit(2);
      }
            
      /* Read next leaf */
      getline(ifs2 ,current_leaf);
      
    } 
  } else {
    cout<<"\nUnable to open leaves-to-be-added file !";
    exit(1);
  }
   
    
  cout<<"\nThe number of leaves to be added is "<<leafcount;
  cout<<"\nThe number of replicates is "<<replicates;  
      
  /* Create replicates of new tree */
  for(int k=0; k<replicates; k++) {    
    int currnodecount = initial_currnodecount;
    int* parent_array = new int[max_total_nodes];
    int* leafindex  = new int[leafcount]; 
    int* added_leaf = new int[leafcount]; 
      
    /* Copy initial parent array and translate array */
    for(int i=0; i<currnodecount; i++) {
      parent_array[i] = initial_parent_array[i];
    }
  
    /* Declare tree, labels & weights */
    aw::Tree t = initial_t;
    aw::idx2name t_name = initial_t_name;
    aw::idx2weight_double t_weight = initial_t_weight;
          
    /* Number of leaves left to be added - temp var */
    int templc = leafcount;
      
    /* Initialise all leaves as not added and copy leaf indices */
    for(int j=0; j<leafcount; j++) {
      leafindex[j] = j;
      added_leaf[j] = NOTADDED;
    }
      
    /* Pick random leaf to add */
    for(int j=0; j<leafcount; j++) {
      int random_leaf;
      int rnum;
      do {
	rnum = rand() % templc;
	random_leaf = leafindex[rnum];      
      } while(added_leaf[random_leaf]==ADDED);
   
      /* Mark current leaf as added */
      added_leaf[random_leaf] = ADDED;     
      
      /* Reduce length of leafindex by 1 and swap chosen leaf with last element */
      leafindex[rnum] = leafindex[templc-1];
      templc--;     
      string newleaf = leaves_array[random_leaf];
      
      /* Get MRCA of current leaves ancestors */
      int lca_id = leaf_lca_array[random_leaf];
      
      /* Find root and parent of root of current subtree */
      unsigned int subtree_root_node = lca_id;
      unsigned int subtree_root_node_parent = parent_array[lca_id];
      
      /* Store branch lengths of current subtree in an array */
      double curr_subtree_bl = 0;
      int subtree_nodecount = 0;
      double* bl_array = new double[max_total_nodes];
      
      /*Traverse subtree of current leaf - included from Andre's email*/
      for (aw::Tree::iterator_postorder v=t.begin_postorder(subtree_root_node,subtree_root_node_parent),vEE=t.end_postorder(); v!=vEE; ++v) {
	unsigned int current_node = v.idx;
	curr_subtree_bl += t_weight[current_node][0];
	translate_index[subtree_nodecount] = current_node;
	bl_array[subtree_nodecount++] = curr_subtree_bl;                
      } 
      
      /* Select a random value and find corresponding edge */
      int untranslated_node;
      unsigned int selected_edge;
      double random_blength;
      do{
        random_blength = (double)rand() * (double)curr_subtree_bl / (double)RAND_MAX; 
        untranslated_node = binarysearch(bl_array, subtree_nodecount, random_blength); 
        selected_edge = translate_index[untranslated_node];         
      }while(selected_edge == subtree_root_node && selected_edge == t.root);
     
      
      /* Obtain individual lengths by subtracting cumulative lengths */
      double original_length = bl_array[untranslated_node] - bl_array[untranslated_node-1];
      random_blength = random_blength - bl_array[untranslated_node-1];
      double reduce_length =   original_length - random_blength; 
            
      /* Change length of branch where inserted */
      t_weight[selected_edge][0] = random_blength;
            
      /* Insert new internal node to attach new leaf */
      unsigned int i_n = t.new_node();
      t_weight[i_n][0] = reduce_length;    
      //bl_array[currnodecount] = bl_array[currnodecount-1] + reduce_length;
  
      /* Insert new leaf */
      unsigned int l_n = t.new_node(); 
      t_weight[l_n][0] = random_blength;
      t_name[l_n] = newleaf; 
      parent_array[l_n] = i_n;
                  
      /* Remove edge b/w selected node and parent */
      unsigned int current_parent = parent_array[selected_edge];
      //cout<<"\nCurrent parent "<<current_parent<<" insert branch "<<selected_edge;
      t.remove_edge(current_parent, selected_edge); 
           
      /* Add edge b/w new leaf and new internal */       
      t.add_edge(i_n, l_n);
          
      /* Add edge b/w new internal and selected node */
      t.add_edge(i_n, selected_edge);
      parent_array[selected_edge] = i_n;
          
      /* Assign parent of new internal to old parent of selected */
      t.add_edge(i_n, current_parent);
      parent_array[i_n] = current_parent;
      delete [] bl_array; 
               
    }
     
    /* Find new total branch length */
    double bl_temp =0;
    TREE_POSTORDER2(node,t){
      unsigned int currnode = node.idx;    
      bl_temp += t_weight[currnode][0]; 
    }   
    cout<<"\nThe new total branch length is: "<<bl_temp;
    
    /*Write tree to file */
    aw::tree2newick(ofs, t, t_name, t_weight);
    ofs<<endl;
   
    delete[] parent_array;
    delete[] leafindex; 
    delete[] added_leaf;
  }

  delete [] initial_parent_array; 
  delete [] translate_index; 
  
  cout<<"\n";
  return 0;
  
  
}  

/* Binary search on larray which is sorted
 * Function returns index of edge to add leaf based on key
 * Returns zero if edge array contains just one element 
*/
int binarysearch(double *larray, int size, double key){
  int first = 0;
  int last = size-1;
  
  /* Single element return 0*/
  if(size == 1) {
    return 0;
  }
  
  /* Binary search */
  while (first <= last) {
    int mid = (first + last) / 2; 
    if (key > larray[mid] && key > larray[mid+1]) 
      first = mid + 1;
    else if (key < larray[mid] && key < larray[mid+1]) 
      last = mid - 1;
    else{
      /* Return next index of edge found*/
      return mid+1;
    }           
  }
  
  /* Insert in zeroth edge */
  return 0;
}

