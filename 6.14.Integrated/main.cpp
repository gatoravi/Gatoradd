/* STEM MRCA Version 1
 * Author - Avinash Ramu using Andre Wehe's TREE Library, Feb 11 2011
 * sample usage -  ./exec tree_file leaves_file replicates opfile
 * This is a program to add taxa of specified family to a tree depending on branch lengths at all edges including stem.
 * usage - ./executable initial_tree_file leaves_families_file no_of_replicates opfile_name 
 * STEM is default, Stem/Crown can be specified by the user
 * Option 1 ( Specify the two species whose LCA is the root of subtree to insert)
 *    TaxonName Species1 Species2 Stem/Crown
 * Option 2 ( Specify the species name)
 *    TaxonName species_name
 * Option 3 ( Specify a part of the Species Name)
 *    TaxonName  Part_of_species_name
 * Option 4 ( Insert anywhere in the tree )
 *    TaxonName RANDOM  
 */
 
 
extern const char *builddate;
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
#include <boost/algorithm/string.hpp> // for stricmp()
#include <math.h>

#define MAX_LEAF_ADD 1000 /* Max Number of new leaves that can be added */ 
#define ADDED 1
#define NOT_ADDED 0
#define MAX_TREE_SIZE 2000 /* Max number of taxa in the input tree */
#define MAX_FAMILY 2000

using namespace std;
int binarysearch(double *larray, int size, double key);

int main(int argc, char* argv[]) { 
  
  if(argc!=5){
    cout<<"Incorrect Arguments ! Exiting! \n\t usage - ./exec tree_file leaves_file replicates opfile\n";
    exit(1);
  }
 
  // set seed of Rand number generator 
  srand(time(NULL));
  
  // ADD OPTION types
  enum{CROWN, STEM, FAM_CROWN, FAM_STEM};
  
  // INVALID LCA ARRAY VALUE
  const int INVALID = -1;
  
  // Read in command line arguments 
  char* treefile = argv[1];
  char* leaves_file = argv[2];  
  unsigned int   replicates = atoi(argv[3]);
  char* opfile = argv[4];
  
  // Map name of the families to family id 
  map<std::string, int> fam2id;
  
  // Read initial starting tree from treefile
  std::ifstream ifs;
  ifs.open (treefile);
  aw::Tree initial_t;
  aw::idx2name initial_t_name;
  aw::idx2weight_double initial_t_weight;   
  if(ifs.good()) {       
    if (!aw::stream2tree(ifs, initial_t, initial_t_name, initial_t_weight)) {
      cout<<"Unable to read tree from file! Exiting!";
      exit(1);
    }    	
  } else {
    cout<<"Unable to open tree file! Exiting!\n";
    exit(1);
  }   
  
  // Open output file 
  ofstream ofs(opfile);
  if(!ofs.good()) {
    cout << "unable to open output file!";
    exit(1);
  } 
  
  // Set the maximum total number of leaves 
  const unsigned int MAX_TOTAL_NODES = 2*(MAX_TREE_SIZE + MAX_LEAF_ADD -1);
  
  //Parent array and translate array for nodes 
  int*    initial_parents = new int[MAX_TOTAL_NODES];
  unsigned int* translate_index = new unsigned int[MAX_TOTAL_NODES];
  
  // Number of leaves and nodes(internal + leaf) 
  int leafn = 0;
  int initial_nodecount = 0;
    
  // Map name of the leaves to their id 
  map<std::string, int> name2id;
  
  // Create a vector of all leaf names 
  std::vector<std::string> leaves; 
  
  double initialBL = 0;
  
  // Traverse the initial tree to create parents array  and leaves vector
  TREE_POSTORDER2(k,initial_t){
    unsigned int currnode = k.idx;
    //cout<<" "<<currnode<<initial_t_name[currnode];
    if(initial_t.is_leaf(currnode)) {
      string leaf = initial_t_name[currnode];
      name2id[leaf] = currnode;
      leaves.push_back(leaf);
      leafn++;
    }    
    initial_nodecount++;
    initial_parents[currnode] = k.parent;
    initialBL += initial_t_weight[currnode][0];
  } 
  
  
  cout<<"\nThe number of leaves in the initial tree is "<<leafn;
  cout<<"\nThe initial total Branch Length is "<<initialBL;
  
  
  // Create an lca object to find lcas of nodes.
  aw::LCA lca;
  lca.create(initial_t);
     
  //Map of family for taxa
  map<std::string, std::string>families;
    
  // Create an array to store the lca values of individual nodes
  int initial_lcas_index[MAX_TOTAL_NODES];
  int initial_lcas[MAX_TOTAL_NODES];
  int ADDoption[MAX_LEAF_ADD];//the option specified for the leaf
  string leaves_array[MAX_LEAF_ADD];//store all leaves to be added
    
  // Read in leaves to be added from file 
  int leafcount = 0;  //number of leaves to be added.
  std::ifstream ifs2(leaves_file);
  string current_leaf;
  if(!ifs2.is_open()) {
    cout<<"\nUnable to open leaves to be added file !";
    exit(1);
  }
  else {
    getline(ifs2 ,current_leaf);
    cout<<"\n\nReading in Leaves and LCAs";    
    while(ifs2.good()) {  
        
      string token1, token2, token3, token4;
      
      // Parse the leaf and options 
      istringstream iss(current_leaf);
      getline(iss, token1, '\t');
      getline(iss, token2, '\t');
      getline(iss, token3, '\t');
      getline(iss, token4, '\t');
      
      cout<<"\n Token 1 "<<token1<<" Token2 "<<token2<<" Token3 "<<token3<<" Token4 "<<token4<<" DONE";
      
      string taxa, base1, base2, familyName;
      
      // ERROR
      if(token2 == "") {
        std::cout<<"\nInvalid Input ! No family or leaves specified for : "<<token1<<" Exiting !\n";
        exit(1);
      }
      
      //------- RANDOM ---------
      else if(token2 == "RANDOM") {
        cout<<"\nRANDOM option";
        ADDoption[leafcount] = CROWN;//anywhere within the tree, so mark it as CROWN
        taxa = token1;
        base1 = initial_t.root;
        base2 = initial_t.root; 
        initial_lcas_index[leafcount] = INVALID;
      } 
      
      //------------ FAMILY ---------------------
      else if(token3 == "" || token3 == "CROWN" || token3 == "STEM") {
        taxa = token1;
        familyName = token2;
        families[taxa] = familyName;
        //select ADDoption, default is STEM
        if(token3 == "CROWN") {
          cout<<"\nFAMILY CROWN option";
          ADDoption[leafcount] = FAM_CROWN;
        }
        else {
          cout<<"\nFAMILY STEM option";
          ADDoption[leafcount] = FAM_STEM;
        }
        initial_lcas_index[leafcount] = INVALID;        
      }
      
      //----------- CROWN -------------------
      else if (token4 == "CROWN" || token4 == "crown") {
       
         taxa = token1;
         base1 = token2;
         base2 = token3; 
         ADDoption[leafcount] = CROWN;
         
         // Get the IDs of the current leaf's bases
         std::map<std::string, int>::iterator index1  = name2id.find(base1);// Note - I use find to check if bases exist
         std::map<std::string, int>::iterator index2  = name2id.find(base2);
         
         //check if the bases exist
         if( index1 == name2id.end() || index2 == name2id.end()) {
           cout<<"\nCannot find bases for "<<taxa<<" in initial tree. Bases are "<<base1<<" , "<<base2<<"\n\tExiting !"<<endl;
           exit(1);
         }
         
         int base1_id = index1->second;
         int base2_id = index2->second;
         
         // Find the MRCA of the ancestors of curr leaf 
         int root = lca.lca(base1_id, base2_id);
         
         // Point the leaf to the index in the lca array 
         initial_lcas_index[leafcount] = root;
         
         // Store the value of lca in the lca array at specified index 
         initial_lcas[root] = root; //This will be modified later if the new leaf is added to the stem of this family
         cout<<"\nCROWN";
      }
      
      // --------------STEM--------------------
      else {
         taxa = token1;
         base1 = token2;
         base2 = token3; 
         ADDoption[leafcount] = STEM;
         
         // Get the IDs of the current leaf's bases
         std::map<std::string, int>::iterator index1  = name2id.find(base1);
         std::map<std::string, int>::iterator index2  = name2id.find(base2);
         
         //check if the bases exist
         if( index1 == name2id.end() || index2 == name2id.end()) {
           cout<<"\nCannot find bases for "<<taxa<<" in initial tree. Bases are "<<base1<<" , "<<base2<<"\n\tExiting !"<<endl;
           exit(1);
         }
         
         int base1_id = index1->second;
         int base2_id = index2->second;
         
         // Find the MRCA of the ancestors of curr leaf 
         int root = lca.lca(base1_id, base2_id);
         
         //Point the leaf to the index in the lca array 
         initial_lcas_index[leafcount] = root;
         
         // Store the value of lca in the lca array at specified index 
         initial_lcas[root] = root; //This will be modified later if the new leaf is added to the stem of this family
         cout<<"\nSTEM";
      }
      
      // Store the leaf label 
      leaves_array[leafcount] = taxa;
      leaves.push_back(taxa);
            
      leafcount++;
      
      // Read the next leaf 
      getline(ifs2 ,current_leaf);
      
    } 
  } 
   
    
  cout<<"\nNumber of leaves to be added is "<<leafcount;
  cout<<"\nNumber of replicates is "<<replicates;
  cout<<"\n";  
     
  // CREATE REPLICATES OF NEW  TREE 
  for(unsigned int k=0; k<replicates; k++) {
  
    cout<<"\n\n\tREPLICATE NUMBER "<<k+1;    
    int nodecount = initial_nodecount;
    
    // Create a parent array and leaf index for each replicate 
    int* parents = new int[MAX_TOTAL_NODES];
    int* leaf_index  = new int[leafcount]; 
    int* added_leaf = new int[leafcount]; 
    
    int lcas_index[MAX_TOTAL_NODES];
    int lcas[MAX_TOTAL_NODES];
      
    // Copy the initial parent array
    for(int i=0; i<nodecount; i++) {
      parents[i] = initial_parents[i];
    }
    
    // copy lca array
    for(int i=0; i<leafcount; i++) {
      if ( ADDoption[i] == STEM || ADDoption[i] == CROWN ) {
          lcas_index[i] = initial_lcas_index[i];
          lcas[initial_lcas_index[i]] =  initial_lcas[initial_lcas_index[i]];
      } 
    }
  
    // Declare tree, labels & weights
    aw::Tree t = initial_t;
    aw::idx2name t_name = initial_t_name;
    aw::idx2weight_double t_weight = initial_t_weight;
          
    // Number of leaves left to be added
    int leaves_remaining = leafcount;
      
    // Initialise all leaves as not added and copy leaf indices
    for(int i=0; i<leafcount; i++) {
      leaf_index[i] = i;
      added_leaf[i] = NOT_ADDED;
    }
      
    // Add Leaves to initial tree 
    for(int j=0; j<leafcount; j++) {
    
      int random_leaf_index;
      int random;
            
      // pick a random LEAF to add to the TREE
      do {
	random = rand() % leaves_remaining;
	random_leaf_index = leaf_index[random];
      } while(added_leaf[random_leaf_index]==ADDED);
   
      added_leaf[random_leaf_index] = ADDED;//mark leaf as added
     
      // shift leaves by 1
      leaf_index[random] = leaf_index[leaves_remaining-1];
      leaves_remaining--;
      string newleaf = leaves_array[random_leaf_index];
      cout<<"\n\n  ADDING "<<newleaf;   
      
      // root of subtree and parent to insert new taxa
      int subtree_root;
      unsigned int subtree_root_parent;
      
      //--------STEM-------
      if( ADDoption[random_leaf_index] == STEM ) {
        cout<<"\nSTEM ADD";         
        int lca_index = lcas_index[random_leaf_index];
        subtree_root = lcas[lca_index];
        subtree_root_parent = parents[subtree_root];
        cout<<"\nTHE STEM CASE ROOT IS "<<subtree_root;
      }
      
      //--------CROWN--------
      else if ( ADDoption[random_leaf_index] == CROWN ) {
        cout<<"\nCROWN ADD";
        int lca_index = lcas_index[random_leaf_index];
        subtree_root = lcas[lca_index];
        subtree_root_parent = parents[subtree_root];
        cout<<"\nTHE CROWN CASE ROOT IS "<<subtree_root;
      }
      
      //--------FAMILY--------
      else if ( ADDoption[random_leaf_index] == FAM_CROWN || ADDoption[random_leaf_index] == FAM_STEM ) {
        
        cout<<"\n\tFAMILY ADD : ";
        string family = families[newleaf];
        cout<<"newleaf is "<<newleaf<<" family is "<<family;
        int flag = 0, left, right;
        
        //LETS FIND LCA HERE      
        TREE_INORDER2(k,t) {
          //cout<<t_name[k.idx];
          string leaf = t_name[k.idx];
          size_t found = leaf.find(family);
          if (found !=string::npos) {
            //cout<<"\nFAMILY MATCH taxon is "<<leaf<<" family is "<<family;            
            if(flag == 0) {
              left = k.idx;
              flag = 1;
            }
            
            else {
              right = k.idx;
            }
            
          }
        }
        
        if(flag == 0) {
          cout<<"\nFamily prefix "<<family<<" not found ! Exiting! \n";
          exit(1);
        }
        
        // Find the MRCA of the ancestors of curr leaf - Create an lca object to find lcas of nodes.
        aw::LCA lca;
        lca.create(t);
        subtree_root = lca.lca(left, right);
        subtree_root_parent = parents[subtree_root];
        cout<<"\n The left base is "<<t_name[left]<<" The right base is "<<t_name[right];
        //cout<<"\n THE FAMILY CASE ROOT IS "<<subtree_root;
        
      }
      
      // Calculate the branch length of the curr family subtree 
      double subtree_bl = 0;
      
      // Store branch lengths of current subtree
      double* bl_array = new double[MAX_TOTAL_NODES];
      int subtree_nodecount = 0;
      
      // Traverse the subtree of current family 
      for (aw::Tree::iterator_postorder v=t.begin_postorder(subtree_root,subtree_root_parent),vEE=t.end_postorder(); v!=vEE; ++v) {
	unsigned int current_node = v.idx;
	//cout<<"\nNode "<<current_node<<" Branch Length: "<<t_weight[current_node][0];
	subtree_bl += t_weight[current_node][0];
	translate_index[subtree_nodecount] = current_node;
	bl_array[subtree_nodecount++] = subtree_bl;
      }
      //cout<<"\n TOTAL SUBTREE BL "<<subtree_bl;
      
      // Select a random branch length and edge
      int untranslated_node;
      int selected_edge; //edge where to add the random leaf
      double randomblength;
      do {
      
	randomblength = (double)rand() * (double)subtree_bl / (double)RAND_MAX; 
	//cout<<"\nRandomblength is "<<randomblength;
	untranslated_node = binarysearch( bl_array, subtree_nodecount, randomblength);
	//cout<<"\nEdge : "<<untranslated_node<<" Translated: "<<translate_index[untranslated_node]; 
	selected_edge = translate_index[untranslated_node];
	
	//Check for root of tree
	if(selected_edge == (signed)t.root) { //STEM option
	  cout<<"\nSelected Root ! continuing !";
	  continue;
	}
	
	//CROWN option check
	else if ( (selected_edge == subtree_root) && ( ADDoption[random_leaf_index] == CROWN || ADDoption[random_leaf_index] == FAM_CROWN) ) { 
	  cout<<"\nSelected subtree root ( not allowed for CROWN ) ! continuing !";
	  continue;
	}
	
	else {
	  break;
	}
	
      } while(1);
            
            
      // Obtain the individual lengths by subtracting the cumulative lengths
      double original_length = bl_array[untranslated_node] - bl_array[untranslated_node-1];
      randomblength = randomblength - bl_array[untranslated_node-1];
      double reduce_length =   original_length - randomblength;       
      
      // Change length of branch where inserted 
      t_weight[selected_edge][0] = randomblength;
            
      // Insert new internal node to attach new leaf 
      unsigned int i_n = t.new_node();
      t_weight[i_n][0] = reduce_length;
  
      // Insert the new leaf 
      unsigned int l_n = t.new_node(); 
      t_name[l_n] = newleaf;
      parents[l_n] = i_n;
      
      // Change LCA ARRAY if selected edge is the root i.e STEM CASE
      if( selected_edge == subtree_root && (ADDoption[random_leaf_index] == STEM || ADDoption[random_leaf_index] == FAM_STEM) ) {
        lcas[subtree_root] = i_n;
      }
      
      // Remove the edge b/w selected node and parent 
      unsigned int current_parent = parents[selected_edge];// parent of selected edge
      t.remove_edge(current_parent, selected_edge);
          
      // Add an edge b/w new internal and selected node 
      t.add_edge(i_n, selected_edge);
      parents[selected_edge] = i_n;
          
      // Assign parent of new internal to old parent of selected 
      t.add_edge(i_n, current_parent);
      parents[i_n] = current_parent;
            
      double leafLength = 0;
      unsigned int parent = parents[i_n];
      
      // Find branch-length of new leaf for ultra-metric tree, DFS from selected edge to leaf
      for (aw::Tree::iterator_dfs v=t.begin_dfs(i_n, parent),vEE=t.end_dfs(); v!=vEE; ++v) {
	unsigned int node = v.idx;
	if (node == i_n)//don't add current nodes length 
	  continue;
	leafLength += t_weight[node][0];
	if(t.is_leaf(node)){
	  break;	  
	}	  
      }        
      
      // Add edge b/w new leaf and new internal 
      t.add_edge(i_n, l_n);
            
      // Weight of new leaf node
      t_weight[l_n][0] = leafLength;       
     
      delete [] bl_array; 
      
      //ostringstream os;
      //aw::tree2newick(os, t, t_name, t_weight);
      //cout<<"\n"<<os.str()<<endl;
         
    }
    // Calculate new branch length   
    double bl_temp =0;
    TREE_POSTORDER2(k,t) {
	unsigned int currnode = k.idx;    
	bl_temp += t_weight[currnode][0]; 
    }   
    cout<<"\nThe new total branch length is: "<<bl_temp;
        
    // Write the tree to file 
    aw::tree2newick(ofs, t, t_name, t_weight);
    ofs<<endl;   
    
    delete[] parents;
    delete[] leaf_index; 
    delete[] added_leaf;
  }
 
  delete [] initial_parents; 
  delete [] translate_index;   
  cout<<"\n";
  return 0;
  
  
}  

// Use Binary search to find the branch where the random value fits in
int binarysearch(double *larray, int size, double key) {

  int first = 0;
  int last = size-1;
  if(size == 1) {
   return 0;
  }
  while (first <= last) {
    int mid = (first + last) / 2;  // compute mid point.
    if (key > larray[mid] && key > larray[mid+1]) 
      first = mid + 1;  // repeat search in top half.
    else if (key < larray[mid] && key < larray[mid+1]) 
      last = mid - 1; // repeat search in bottom half.
    else{
      return mid+1; // found it. return position 
    }           
  }
  return 0;
  
}
