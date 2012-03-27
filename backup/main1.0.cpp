/* STEM MRCA Version 1
 * Written by Avinash using Andre Wehe's Library Feb 11 2011
 * sample usage -  ./exec tree_file leaves_file replicates opfile
 * This is a program to add taxa of specified family to a tree .
 * usage - ./executable initial_tree_file leaves_families_file no_of_replicates opfile_name 
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
#include <math.h>

#define MAXADDLEAF 1000 /* Max Number of new leaves that can be added */ 
#define ADDED 1
#define NOTADDED 0
#define MAXTREETAXA 2000 /* Max number of taxa in the input tree */

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
  }
  else {
    cout<<"Unable to open tree file! Exiting!\n";
    exit(1);
  }   
  
  /* Open output file */
  ofstream ofs(opfile);
  if(!ofs.good()) {
    cout<<"unable to open output file!";
    exit(1);
  } 
  
  /* Print out tree */
  //ostringstream os;
  //aw::tree2newick(os, t, t_name, t_weight);
  //cout<<"\nOutput : "<<os.str()<<endl;
  
  /* Create a set of all families */
  set<string> family_set;
  
  
  int max_total_nodes = 2*(MAXTREETAXA + MAXADDLEAF -1);
  
  int*    initial_parent_array = new int[max_total_nodes];
  unsigned int* translate_index = new unsigned int[max_total_nodes];
  //double total_blength = 0;
  int leafn = 0;
  int initial_currnodecount = 0;
    
  /* Map the name of the leaves to their id */
  map<std::string, int> name2id_map;
  
  /* Traverse the tree to create parents array  */
  cout<<"\n POST ORDER ";
  TREE_POSTORDER2(k,initial_t){
    unsigned int currnode = k.idx;
    cout<<"\nNode "<<currnode<<" "<<initial_t_name[currnode];
    //cout<<" is "<<t_weight[currnode][0];    
    //total_blength += t_weight[currnode][0]; 
    //translate_index[currnodecount] = currnode;
    //if(currnode != t.root)
      //bl_array[currnodecount] = total_blength;//this array stores cumulative branch lengths & hence is sorted.
    if(initial_t.is_leaf(currnode)) {
        name2id_map[initial_t_name[currnode]] = currnode;
        leafn++;
        //cout<<endl<<t_name[currnode];
    }    
    initial_currnodecount++;
    initial_parent_array[currnode] = k.parent;     
  } 
  cout<<"\nThe number of leaves in the initial tree is "<<leafn;
  cout<<"\nThe number of internal nodes in the initial tree is "<<initial_currnodecount;
  
  /*
  show content of the map 
  map<std::string, int>::iterator it;
  string temp;
  for ( it=name2id_map.begin() ; it != name2id_map.end(); it++ ) {
    cout << (*it).first << " => " << (*it).second << endl;
    temp = (*it).first;
  }
  */
 
  
  /*cout<<"\n\nBEFORE ADDITION\n\n";
  for(int i=0; i<currnodecount-1; i++){
    int translate = translate_index[i]; 
    cout<<"\nbranch number "<<translate;
    if(i!=0)cout<<" length "<<bl_array[i] - bl_array[i-1];
    else cout<<" length "<<bl_array[i];
    cout<<" cumul length "<<bl_array[i];
  }*/
  
  /* Create an lca object to find lcas of nodes. */
  aw::LCA lca;
  lca.create(initial_t);
  
  /* Create an array to store the lca values of induv nodes */
  int leaf_lca_array[MAXADDLEAF];
  bool single_taxa_family[MAXADDLEAF];
  string leaves_array[MAXADDLEAF];
  
  /* Read in leaves to be added from file */
  int leafcount = 0;  //number of leaves to be added.
  std::ifstream ifs2(leaves_file);
  string current_leaf;
  int same = 0;
  if(ifs2.is_open()){
    getline(ifs2 ,current_leaf);
    string taxa, anc1, anc2;
    pair<set<string>::iterator, bool> family_iterator;
    
    while(ifs2.good()) {      
    
      /* Parse the leaf to be added and its two ancestors */
      istringstream iss(current_leaf);
      getline(iss, taxa, '\t');
      getline(iss, anc1, '\t');
      getline(iss, anc2, '\t');
      //cout<<endl<<endl<<current_leaf<<endl;
      //cout<<"Taxa "<<taxa<<endl<<"anc1 is "<<anc1<<endl<<" anc2 is "<<anc2<<endl; 
      
      /* Add the family of current leaf to the family set */
      string curr_leaf_family = anc1 + " " + anc2;
      family_iterator = family_set.insert(curr_leaf_family);
           
      /* Get the IDs of the current leaf's ancestors*/
      int anc1_id = name2id_map[anc1];
      int anc2_id = name2id_map[anc2];
      //cout<<endl<<"anc1_id is: "<<anc1_id;
      //cout<<endl<<"anc2_id is: "<<anc2_id;
      
      /*Identify families with single taxa */
      if(anc1_id == anc2_id) {
        single_taxa_family[leafcount] = true;
        same++;
      } 
      else 
        single_taxa_family[leafcount] = false;
      
      /* Find the MRCA of the ancestors of curr leaf */
      int lca_id = lca.lca(anc1_id, anc2_id);
      //cout<<endl<<"The MRCA of curr leaf anc is: "<<lca_id;
      
      /* Store the leaf label */
      leaves_array[leafcount] = taxa;
      
      /* Store the lca in the array */
      leaf_lca_array[leafcount++] = lca_id;
      
      /* Read the next leaf */
      getline(ifs2 ,current_leaf);
      
    } 
   } 
   else {
    cout<<"\nUnable to open leaves to be added file !";
    exit(1);
   }
   
    
    cout<<"\nThe number of leaves to be added is "<<leafcount;
    cout<<"\nThe number of replicates is "<<replicates;  
    cout<<"\nInitial number of nodes is "<<initial_currnodecount;
      
    /* CREATE REPLICATES OF NEW  TREE */
    for(int k=0; k<replicates; k++) {
    
      int currnodecount = initial_currnodecount;
      //cout<<"Initial curr node count is: "<<currnodecount;
      //double total_blength = total_blength_initial;
      cout<<"\n\nReplicate number : "<<k+1;
      //cout<<"\nThe total branch length is "<<total_blength;
      //double* bl_array = new double[total_nodes];
      int* parent_array = new int[max_total_nodes];
      int* leafindex  = new int[leafcount]; 
      int* added_leaf = new int[leafcount]; 
      //unsigned int* translate_index = new unsigned int[total_nodes];
      
      /* Copy the initial parent array and translate array */
      for(int i=0; i<currnodecount; i++) {
        //translate_index[i] = translate_index_initial[i];
        parent_array[i] = initial_parent_array[i];
        //bl_array[i] = bl_array_initial[i];
      }
  
      /* Declare tree, labels & weights */
      aw::Tree t = initial_t;
      aw::idx2name t_name = initial_t_name;
      aw::idx2weight_double t_weight = initial_t_weight;
    
      //aw::Tree t;
      //aw::idx2name t_name;
      //aw::idx2weight_double t_weight;
      //aw::stream2tree(ifs, t, t_name, t_weight);
      
      /* Number of leaves left to be added */
      int templc = leafcount;
      
      /* Initialise all leaves as not added and copy leaf indices */
      for(int j=0; j<leafcount; j++) {
        leafindex[j] = j;
        added_leaf[j] = NOTADDED;
      }
      
      for(int j=0; j<leafcount; j++) {
      int random_leaf;
      int rnum;
      do {
	rnum = rand() % templc;
	random_leaf = leafindex[rnum];      
	//cout<<"\nj"<<j;
	//cout<<"\nTEMPLC"<<templc;
	//cout<<"\nleafcount"<<leafcount;
	//cout<<endl<<"rnum "<<rnum<<"random "<<random_leaf<<endl<<" "<<"added y/n "<<added_leaf[random_leaf];
	//cout<<"chk";
      }while(added_leaf[random_leaf]==ADDED);
   
      added_leaf[random_leaf] = ADDED;
      //cout<<"\nLeaf Number "<<random_leaf+1; // leaf numbers 1 to n
     
      /* Reduce length of leafindex by 1 and swap chosen leaf with last element */
      leafindex[rnum] = leafindex[templc-1];
      templc--;     
      string newleaf = leaves_array[random_leaf];
      
      /* Get the MRCA of the current leaves ancestors */
      int lca_id = leaf_lca_array[random_leaf];
      
      unsigned int subtree_root_node = lca_id;
      unsigned int subtree_root_node_parent = parent_array[lca_id];
      cout<<endl<<"\nCurrent leaf family subtree root is "<<subtree_root_node;
      cout<<"\nCurrent leaf family subtree root parent is "<<subtree_root_node_parent<<endl;
      double curr_subtree_bl = 0;
      cout<<"\nTraversal for current leaf's subtree leafnumber: "<<leafcount<<" leafname "<<newleaf;
      int subtree_nodecount = 0;
      double* bl_array = new double[max_total_nodes];
      /*Traverse the subtree of current leaf */
      for (aw::Tree::iterator_postorder v=t.begin_postorder(subtree_root_node,subtree_root_node_parent),vEE=t.end_postorder(); v!=vEE; ++v) {
                unsigned int current_node = v.idx;
                cout<<"\nNode "<<current_node<<" Branch Length: "<<t_weight[current_node][0];
                curr_subtree_bl += t_weight[current_node][0];
                translate_index[subtree_nodecount] = current_node;
                bl_array[subtree_nodecount++] = curr_subtree_bl;                
      }          
      
      cout<<"\ncurr_subtree_bl = "<<curr_subtree_bl;
      int untranslated_node;
      double randomblength;
      do{
        randomblength = (double)rand() * (double)curr_subtree_bl / (double)RAND_MAX; 
        cout<<"\nRandomblength is "<<randomblength;
        untranslated_node = binarysearch(bl_array, subtree_nodecount, randomblength);
        cout<<"\nEdge : "<<untranslated_node<<" Translated: "<<translate_index[untranslated_node];        
      }while(translate_index[untranslated_node] == subtree_root_node );
      
      
      for(int i=0; i<subtree_nodecount; i++){
        cout<<"\ncumul branch length "<<i<<" "<<bl_array[i];
        if(i!=0)cout<<" actual length "<<bl_array[i] - bl_array[i-1];
     }
      
      int add_edge = translate_index[untranslated_node];
      cout<<"\nedge where to add "<<add_edge;
      cout<<"\nuntranslated_node "<<untranslated_node;
      
      double original_length;
      
      /* obtain the individual lengths by subtracting the cumulative lengths */
      original_length = bl_array[untranslated_node] - bl_array[untranslated_node-1];
      randomblength = randomblength - bl_array[untranslated_node-1];
      double reduce_length =   original_length - randomblength; 
      cout<<"\nOriginal Length "<<original_length;
      cout<<"\nrandomblength "<<randomblength;
      cout<<"\nreduce_length "<<reduce_length;
      
      
      /* change length of branch where inserted */
      t_weight[add_edge][0] = randomblength;
      
      
      /* insert new internal node to attach new leaf */
      unsigned int i_n = t.new_node();
      t_weight[i_n][0] = reduce_length;    
      //bl_array[currnodecount] = bl_array[currnodecount-1] + reduce_length;
  
      /* insert the new leaf */
      unsigned int l_n = t.new_node(); 
      t_weight[l_n][0] = randomblength;
      t_name[l_n] = newleaf; 
      parent_array[l_n] = i_n;
      
      
       /* remove the edge b/w selected node and parent */
      unsigned int current_parent = parent_array[add_edge];
      //cout<<"\nCurrent parent "<<current_parent<<" insert branch "<<add_edge;
      t.remove_edge(current_parent, add_edge); 
      
      
      //cout<<"\nnot root selected";
      //cout<<"\nThe parent is : "<<parent;	      
     
      /* add an edge b/w new leaf and new internal */       
      t.add_edge(i_n, l_n);
      
    
      /* add an edge b/w new internal and selected node */
      t.add_edge(i_n, add_edge);
      parent_array[add_edge] = i_n;
      
    
      /* assign parent of new internal to old parent of selected */
      t.add_edge(i_n, current_parent);
      parent_array[i_n] = current_parent;
       
      unsigned int edgec = t.edge_size();
      unsigned int  nodec = t.node_size();
      
      cout<<"\nCurrent number of edges is "<<edgec;
      cout<<"\nCurrent number of nodes is "<<nodec;
      
      /* Output tree to verify */
      TREE_POSTORDER2(k,t){
	unsigned int currnode = k.idx;
	cout<<"\nWeight for node "<<currnode<<" "<<t_name[currnode];
	cout<<" is "<< t_weight[currnode][0];    
	//total_blength += t_weight[currnode][0]; 
	//translate_index[currnodecount] = currnode;
	//if(currnode!=t.root)
	  //bl_array[currnodecount++] = total_blength;//this array stores cumulative branch lengths & hence is sorted.
	//if(t.is_leaf(k.idx)){
        //leafn++;
	//}
	//parent_array[currnode] = k.parent;     
      }   
      //int wait;
      //cin>>wait;
      delete [] bl_array; 
      
         
    }
   ostringstream os;
   aw::tree2newick(os, t, t_name, t_weight);
   cout<<"\nOutput : "<<os.str()<<endl;
   
   aw::tree2newick(ofs, t, t_name, t_weight);
   ofs<<endl;
   //cout<<"\nThe new total branch length is: "<<total_blength;
   
   delete[] parent_array;
   delete[] leafindex; 
   delete[] added_leaf;
  }
  
  //cout<<"\nThe size of the name2id map is "<<name2id_map.size();
  //cout<<"\nThe total number of distinct families is "<<(int) family_set.size();
  //cout<<"\nThe total number of leaves to be added is "<<leafcount;
  //cout<<"\nThe number of families with one taxa is "<<same;
  //int currleafcount = leafn;  
  //int* added_leaf = new int[leafn + leafcount]; 
  //int* leafindex  = new int[leafn + leafcount]; 
  
  
  
  
  
 
  
  /*for (aw::Tree::iterator_postorder v=t.begin_postorder(subtree_root_node,subtree_root_node_parent),vEE=t.end_postorder(); v!=vEE; ++v) {
  
                unsigned int current_node = v.idx;
                cout<<"\nThe current node is "<<current_node;
  }*/
  
  //aw::tree2newick(os, t, t_name, t_weight);
  //cout<<"\nOutput : "<<os.str()<<endl;
  
  
  /* Print out tree */
 
  //cout<<"\nOutput : "<<os.str()<<endl;
  /* deallocate memory */
  //delete [] added_leaf; 
  //delete [] leafindex; 
 
  delete [] initial_parent_array; 
  delete [] translate_index; 
  
  cout<<"\n";
  return 0;
  
  
}  


int binarysearch(double *larray, int size, double key){
  int first = 0;
  int last = size-1;
  //cout<<"\nSize of array "<<size<<" key "<<key<<" last "<<last;
  while (first <= last) {
    int mid = (first + last) / 2;  // compute mid point.
    //cout<<"\nlarraymid "<<mid<<" "<<larray[mid]<<"\nlarraymid+1 "<<mid+1<<" "<<larray[mid+1];
    if (key > larray[mid] && key > larray[mid+1]) 
      first = mid + 1;  // repeat search in top half.
    else if (key < larray[mid] && key < larray[mid+1]) 
      last = mid - 1; // repeat search in bottom half.
    else{
      //cout<<"\nReturned value: "<<mid+1;
      return mid+1;     // found it. return position 
    }
           
  }
  return 0;
}
