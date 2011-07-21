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
//#define FAMILY 0
//#define CROWN 1
//#define STEM 2

using namespace std;
int binarysearch(double *larray, int size, double key);

int main(int argc, char* argv[]) { 
  
  if(argc!=5){
    cout<<"Incorrect Arguments ! Exiting! \n\t usage - ./exec tree_file leaves_file replicates opfile\n";
    exit(1);
  }
 
  /* set seed of Rand number generator */
  srand(time(NULL));
  
  //ADD OPTION types
  enum{CROWN, STEM, FAM_CROWN, FAM_STEM};
  
  /* Read in command line arguments */
  char* treefile = argv[1];
  char* leaves_file = argv[2];  
  unsigned int   replicates = atoi(argv[3]);
  char* opfile = argv[4];
  //unsigned int family_ends[MAX_FAMILY][2];
  
  /* Map name of the families to family id */
  map<std::string, int> fam2id;
  
  /* Read initial starting tree from treefile */
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
  
  /* Open output file */
  ofstream ofs(opfile);
  if(!ofs.good()) {
    cout << "unable to open output file!";
    exit(1);
  } 
  
  /* Print out tree */
  //ostringstream os;
  //aw::tree2newick(os, t, t_name, t_weight);
  //cout<<"\nOutput : "<<os.str()<<endl;
  
  /* Create a set of all species families */
  //set<string> families;
  
  /* Set the maximum total number of leaves */
  unsigned int max_total_nodes = 2*(MAX_TREE_SIZE + MAX_LEAF_ADD -1);
  
  /* Define Parent array and translate array for nodes */
  int*    initial_parent = new int[max_total_nodes];
  unsigned int* translate_index = new unsigned int[max_total_nodes];
  //double total_blength = 0;
  
  /* leaf count and node count */
  int leafn = 0;
  int initial_nodecount = 0;
    
  /* Map name of the leaves to their id */
  map<std::string, int> name2id;
  
  /*Create a vector of all leaf names */
  std::vector<std::string> leaves; 
  
  //unsigned int familyc =0;//number of families
  
  /* Traverse the initial tree to create parents array  and leaves vector*/
  TREE_POSTORDER2(k,initial_t){
    unsigned int currnode = k.idx;
    cout<<" "<<currnode<<initial_t_name[currnode];
    //string family;//family of current leaf
    
    //cout<<"\nNode "<<currnode<<" "<<initial_t_name[currnode];
    //cout<<" is "<<t_weight[currnode][0];    
    //total_blength += t_weight[currnode][0]; 
    //translate_index[currnodecount] = currnode;
    //if(currnode != t.root)
    //bl_array[currnodecount] = total_blength;//this array stores cumulative branch lengths & hence is sorted.
    if(initial_t.is_leaf(currnode)) {
      string leaf = initial_t_name[currnode];
      name2id[leaf] = currnode;
      leaves.push_back(leaf);
      leafn++;
      
      /* Add the family of current leaf to the family set */
      /*istringstream is(leaf);
      getline(is, family, '_');//extract the family name
      pair<set<string>::iterator, bool> family_iterator = families.insert(family);      
      if(family_iterator.second) {//new family added to set
        cout<<"\nThe new family is "<<family<<" the family number is "<<familyc<<" the curr node "<<currnode;
        fam2id[family] = familyc;
        family_ends[familyc][0] = family_ends[familyc][1] = currnode;//the left and right ends of new family
        familyc++;       
      } else {//already added family
        int famNum = fam2id[family];//get the family number
        cout<<"\nThe old family is "<<family<<" the family number is "<<famNum<<" the curr node "<<currnode;
        if ( currnode > family_ends[famNum][1]) {
          family_ends[famNum][1] = currnode;//modify right end
        }
        else if ( currnode < family_ends[famNum][0]) {
          family_ends[famNum][0] = currnode;//modify left end
        }
      }*/     
      //cout<<endl<<t_name[currnode];
      
    }    
    initial_nodecount++;
    initial_parent[currnode] = k.parent;     
  } 
  
  /* Create an lca object to find lcas of nodes. */
  aw::LCA lca;
  lca.create(initial_t);
  
  /*Roots of families */
  //int famroots[MAX_FAMILY];
  
  //cout<<"\nThe total number of families is "<<familyc<<"\n The family ends are ";
  
  /*for ( unsigned int j = 0; j< familyc; j++ ) {
    cout<<"\n Family "<<j<<" "<<family_ends[j][0]<<" "<<family_ends[j][1];
    //Find the MRCA of the ancestors of curr leaf 
    int root = lca.lca(family_ends[j][0], family_ends[j][1]);
    famroots[j] = root;
    cout<<"\n The root of family "<<j<<" is "<<root;
  }*/
  
  cout<<"\nThe number of leaves in the initial tree is "<<leafn;
  //cout<<"\nThe number of internal nodes in the initial tree is "<<initial_nodecount;
  //exit(1);
  /*
    //show content of the map 
   
    string temp;
    for ( it=name2id.begin() ; it != name2id.end(); it++ ) {
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
  
  //Map of left and right base taxa for adding each taxa
  //map<std::string, std::string>left;
  //map<std::string, std::string>left;
  
  //Map of family for taxa
  map<std::string, std::string>families;
  
  
  /* Create an array to store the lca values of individual nodes */
  int initial_leaf_lca_array_index[MAX_LEAF_ADD];
  int initial_lca_array[MAX_LEAF_ADD];
  int ADDoption[MAX_LEAF_ADD];//the option specified for the leaf
  //bool single_taxa_family[MAX_LEAF_ADD];
  string leaves_array[MAX_LEAF_ADD];//store all leaves to be added
    
  /* Read in leaves to be added from file */
  int leafcount = 0;  //number of leaves to be added.
  std::ifstream ifs2(leaves_file);
  string current_leaf;
  //int same = 0;// Number of taxa with single family
  if(!ifs2.is_open()) {
    cout<<"\nUnable to open leaves to be added file !";
    exit(1);
  }
  else {
    getline(ifs2 ,current_leaf);
    
    //pair<set<string>::iterator, bool> family_iterator;
    
    cout<<"\n\nReading in Leaves and LCAs";    
    while(ifs2.good()) {  
        
      string token1, token2, token3, token4;
      
      /* Parse the leaf and options */
      istringstream iss(current_leaf);
      getline(iss, token1, '\t');
      getline(iss, token2, '\t');
      getline(iss, token3, '\t');
      getline(iss, token4, '\t');
      
      cout<<"\n Token 1 "<<token1<<" Token2 "<<token2<<" Token3 "<<token3<<" Token4 "<<token4<<" DONE";
      //cout<<endl<<endl<<current_leaf<<endl;
      //cout<<"Taxa "<<taxa<<endl<<"base1 is "<<base1<<endl<<" base2 is "<<base2<<endl; 
      
      string taxa, base1, base2, familyName;
      
      if(token2 == "") {//ERROR
        std::cout<<"\nInvalid Input ! No family or leaves specified for : "<<token1<<" Exiting !\n";
        exit(1);
      }
      //RANDOM
      else if(token2 == "RANDOM") {
        cout<<"\nRANDOM option";
        ADDoption[leafcount] = CROWN;//anywhere within the tree, so mark it as CROWN
        taxa = token1;
        base1 = initial_t.root;
        base2 = initial_t.root; 
      } 
      
      //FAMILY
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
        
        //int familyLen = familyName.length();
        
        //ADDoption[leafcount] = FAMILY;
        //cout<<"\nfamily ADD";
        
        /*int max = 0;
        int min = numeric_limits<int>::max();
        TREE_INORDER2(k,t){
          cout<<t_name[k.idx];
          string leaf = t_name[k.idx];
          size_t found = leaf.find(family);
          if (found !=string::npos) {
          
            cout<<"\n FAMILY MATCH taxon is "<<leaves[i];
            
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
          cout<<"FAMILY "<<family<<" NOT FOUND! EXITING ! \n";
          exit(1);
        }*/
        //Find the two end taxa for this leaf-family 
        /*for(unsigned int i=0; i<leaves.size(); i++) {
          //string leaf = leaves[i];
          size_t found = leaves[i].find(familyName);
          if (found !=string::npos) {
            cout<<"\n FAMILY MATCH taxon is "<<leaves[i];
            
            if(name2id[leaves[i]] < min ) {
              base1 = leaves[i];
              min = name2id[leaves[i]];//update min
            }
            
            if(name2id[leaves[i]] > max ) {
              base2 = leaves[i];
              max = name2id[leaves[i]];//update max
            }
            
          }
        }*/
        
        
        
        
        //cout<<"\nANC1 "<<base1<<" base2 "<<base2;
        
        //cout<<"\n Exiting !";
        //exit(1);
        
        /*int num = fam2id[familyName];// 7/15/2011 
        cout<<"\n The family number is "<<num;
        cout<<"\n The root of the family is "<<famroots[num];*/
        
        
        
        /*//Find family for leaf
        map<std::string, int>::iterator it;
        for ( it=name2id.begin() ; it != name2id.end(); it++ ) {
          cout << "\n" << (*it).first << " => " << (*it).second;
          string leaf = (*it).first;
          size_t found;
          found=leaf.find(familyName);
          if (found!=string::npos) 
            cout<<"\nFOUND FAMILY LEAF "<<leaf;
        }*/        
      
        
      }
      //CROWN
      else if (token4 == "CROWN" || token4 == "crown") { 
         taxa = token1;
         base1 = token2;
         base2 = token3; 
         ADDoption[leafcount] = CROWN;
         /* Get the IDs of the current leaf's ancestors*/
         int base1_id = name2id[base1];
         int base2_id = name2id[base2];
         /* Find the MRCA of the ancestors of curr leaf */
         int root = lca.lca(base1_id, base2_id);
         //Point the leaf to the index in the lca array */
         initial_leaf_lca_array_index[leafcount] = root;
         /* Store the value of lca in the lca array at specified index */
         initial_lca_array[root] = root; //This will be modified later if the new leaf is added to the stem of this family
         cout<<"\nCROWN";
      }
      // STEM
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
           cout<<"\nCannot find bases for "<<taxa<<" in initial tree. Bases are "<<token2<<" , "<<token3<<"\n\tExiting !"<<endl;
           exit(1);
         }
         
         int base1_id = index1->second;
         int base2_id = index2->second;
         
         /* Find the MRCA of the ancestors of curr leaf */
         int root = lca.lca(base1_id, base2_id);
         
         //Point the leaf to the index in the lca array */
         initial_leaf_lca_array_index[leafcount] = root;
         
         /* Store the value of lca in the lca array at specified index */
         initial_lca_array[root] = root; //This will be modified later if the new leaf is added to the stem of this family
         cout<<"\nSTEM";
      }
      
      /* Store the leaf label */
      leaves_array[leafcount] = taxa;
      leaves.push_back(taxa);
      
      
      //if(family_iterator.second) {
        //cout<<"\n NEW"<<leafcount;
      //} else {
        //cout<<"\nOld !!"<<leafcount;
      //}
      
      
      /*Identify families with single taxa */
      
      
      /* Find the MRCA of the ancestors of curr leaf */
      //int root = lca.lca(base1_id, base2_id);
      //cout<<endl<<"The MRCA of curr leaf anc is: "<<lca_id;
            
      /* Point the leaf to the index in the lca array */
      //initial_leaf_lca_array_index[leafcount++] = root;
      
      /* Store the value of lca in the lca array at specified index */
      //initial_lca_array[root] = root; //This will be modified later if the new leaf is added to the stem of this family
           
      /*cout<<"\nLeaf number "<<leafcount+1<<" base1_id is: "<<base1_id<<" base2_id is: "<<base2_id<<" lca "<<root;
      if(base1_id == base2_id) {
        cout<<" same 2 Ancestors.";
        //exit(1);
      } 
      if(root == 0) {
        cout<<" LCA is root !\n";
        //exit(2);
      }*/
      
      
      /*Identify families with single taxa */
      //if(base1_id == base2_id) {
        //single_taxa_family[leafcount] = true;
        //cout<<"\n STF base1 "<<base1_id<<" base2 "<<base2_id<<" lca "<<lca_id;
        //same++;
      //} 
      //else 
        //single_taxa_family[leafcount] = false;
      
      
      
      leafcount++;
      
      /* Read the next leaf */
      getline(ifs2 ,current_leaf);
      
    } 
  } 
   
    
  cout<<"\nNumber of leaves to be added is "<<leafcount;
  cout<<"\nNumber of replicates is "<<replicates;  
  //cout<<"\nInitial number of nodes is "<<initial_nodecount;
  cout<<"\n";  
  
  //exit(1);
  
  
  
   
  /* CREATE REPLICATES OF NEW  TREE */
  for(unsigned int k=0; k<replicates; k++) {
    cout<<"\n\n\tREPLICATE NUMBER "<<k+1;
    //int wait;
    //cin>>wait;
    
    int currnodecount = initial_nodecount;
    //cout<<"Initial curr node count is: "<<currnodecount;
    //double total_blength = total_blength_initial;
    //cout<<"\n\nReplicate number : "<<k+1;
    //cout<<"\nThe total branch length is "<<total_blength;
    //double* bl_array = new double[total_nodes];
    
    /*Create a parent array and leaf index for each replicate */
    int* parent_array = new int[max_total_nodes];
    int* leafindex  = new int[leafcount]; 
    int* added_leaf = new int[leafcount]; 
    int leaf_lca_array_index[MAX_LEAF_ADD];
    int lca_array[MAX_LEAF_ADD];
      
    /* Copy the initial parent array  */
    for(int i=0; i<currnodecount; i++) {
      parent_array[i] = initial_parent[i];
      leaf_lca_array_index[i] = initial_leaf_lca_array_index[i];
      lca_array[i] =  initial_lca_array[i];
    }
  
    /* Declare tree, labels & weights */
    aw::Tree t = initial_t;
    aw::idx2name t_name = initial_t_name;
    aw::idx2weight_double t_weight = initial_t_weight;
          
    /* Number of leaves left to be added */
    int templc = leafcount;
      
    /* Initialise all leaves as not added and copy leaf indices */
    for(int j=0; j<leafcount; j++) {
      leafindex[j] = j;
      added_leaf[j] = NOT_ADDED;
    }
      
    /* Add Leaves to initial tree */
    for(int j=0; j<leafcount; j++) {
    
      //cout<<"\nAdding leafnumber "<<j+1;
      //int wait;
      //cin>>wait;
           
      int random_leaf;
      int rnum;
      
      
      /* Pick a random LEAF to add to the TREE*/
      do {
	rnum = rand() % templc;
	random_leaf = leafindex[rnum];
	//cout<<"\nj"<<j;
	//cout<<"\nTEMPLC"<<templc;
	//cout<<"\nleafcount"<<leafcount;
	//cout<<endl<<"rnum "<<rnum<<"random "<<random_leaf<<endl<<"  "<<"added y/n "<<added_leaf[random_leaf];
	//cout<<"chk";
      } while(added_leaf[random_leaf]==ADDED);
   
      added_leaf[random_leaf] = ADDED;//mark leaf as added
      //cout<<"\nLeaf Number "<<random_leaf+1; // leaf numbers 1 to n
     
      // shift leaves by 1
      leafindex[rnum] = leafindex[templc-1];
      templc--;
      string newleaf = leaves_array[random_leaf];   
      cout<<"\n\n  ADDING "<<newleaf;   
      
      //root of subtree to insert new taxa
      int subtree_root_node;
      unsigned int subtree_root_node_parent;
      
      //STEM
      if( ADDoption[random_leaf] == STEM ) {
        cout<<"\nSTEM ADD";         
        int lca_id_index = leaf_lca_array_index[random_leaf];
        subtree_root_node = lca_array[lca_id_index];
        subtree_root_node_parent = parent_array[subtree_root_node];
        cout<<"\nTHE STEM CASE ROOT IS "<<subtree_root_node;
      }
      
      //CROWN
      else if ( ADDoption[random_leaf] == CROWN ) {
        cout<<"\nCROWN ADD";
        int lca_id_index = leaf_lca_array_index[random_leaf];
        subtree_root_node = lca_array[lca_id_index];
        subtree_root_node_parent = parent_array[subtree_root_node];
        cout<<"\nTHE CROWN CASE ROOT IS "<<subtree_root_node;
      }
      
      //FAMILY
      else if ( ADDoption[random_leaf] == FAM_CROWN || ADDoption[random_leaf] == FAM_STEM ) {
        
        cout<<"\nFAMILY ADD";
        string family = families[newleaf];
        cout<<"\n newleaf is "<<newleaf<<" family is "<<family;
        int flag = 0, left, right;
        
        //LETS FIND LCA HERE      
        TREE_INORDER2(k,t){
          //cout<<t_name[k.idx];
          string leaf = t_name[k.idx];
          size_t found = leaf.find(family);
          if (found !=string::npos) {          
            cout<<"\nFAMILY MATCH taxon is "<<leaf<<" family is "<<family;            
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
        subtree_root_node = lca.lca(left, right);
        subtree_root_node_parent = parent_array[subtree_root_node];
        cout<<"\n The left base is "<<t_name[left]<<" The right base is "<<t_name[right];
        cout<<"\n THE FAMILY CASE ROOT IS "<<subtree_root_node;
        
      }
      
      
      
      
      
      
      ///* Get the MRCA of the current leaves ancestors */
      //int lca_id_index = leaf_lca_array_index[random_leaf];
      //int lca_id = lca_array[lca_id_index];
      
      /* The root of the family subtree is the lca of the two leaves in the ip file */
      //unsigned int subtree_root_node = lca_id;
     
      //cout<<endl<<"\nCurrent leaf "<<j+1<<" family subtree root is "<<subtree_root_node<<" root parent is "<<subtree_root_node_parent<<endl;
      
      /* Calculate the branch length of the curr family subtree */
      double curr_subtree_bl = 0;
      //cout<<"\nTraversal for current leaf's subtree leafnumber: "<<leafcount<<" leafname "<<newleaf;
      
      // Store branch lengths of current subtree
      double* bl_array = new double[max_total_nodes];
      int subtree_nodecount = 0;
      
      //Traverse the subtree of current family 
      for (aw::Tree::iterator_postorder v=t.begin_postorder(subtree_root_node,subtree_root_node_parent),vEE=t.end_postorder(); v!=vEE; ++v) {
	unsigned int current_node = v.idx;
	//cout<<"\nNode "<<current_node<<" Branch Length: "<<t_weight[current_node][0];
	curr_subtree_bl += t_weight[current_node][0];
	translate_index[subtree_nodecount] = current_node;
	bl_array[subtree_nodecount++] = curr_subtree_bl;
      }
      
      // Select a random branch length and edge
      //cout<<"\ncurr_subtree_bl = "<<curr_subtree_bl;
      int untranslated_node;
      int selected_edge; //edge where to add the random leaf
      double randomblength;
      do {
	randomblength = (double)rand() * (double)curr_subtree_bl / (double)RAND_MAX; 
	cout<<"\nRandomblength is "<<randomblength;
	untranslated_node = binarysearch( bl_array, subtree_nodecount, randomblength);
	cout<<"\nEdge : "<<untranslated_node<<" Translated: "<<translate_index[untranslated_node]; 
	selected_edge = translate_index[untranslated_node];
	
	//Check for root of tree/subtree
	if(selected_edge == (signed)t.root) { //STEM option
	  cout<<"\nSelected Root ! continuing !";
	  continue;
	}
	
	//CROWN option check
	else if ( (selected_edge == subtree_root_node) && ( ADDoption[random_leaf] == CROWN || ADDoption[random_leaf] == FAM_CROWN) ) { 
	  cout<<"\nSelected subtree root ( not allowed for CROWN ) ! continuing !";
	  continue;
	}
	else {
	  break;
	}
	
      } while(1);
            
      /* Print cumulative branch lengths and individual branch lengths */
      /*for(int i=0; i<subtree_nodecount; i++){
        cout<<"\ncumul branch length "<<i<<" "<<bl_array[i];
        if(i!=0)cout<<" actual length "<<bl_array[i] - bl_array[i-1];
      }
      */    
      //cout<<"\nuntranslated_node "<<untranslated_node<<" edge where to add "<<selected_edge;
            
      /* obtain the individual lengths by subtracting the cumulative lengths */
      double original_length = bl_array[untranslated_node] - bl_array[untranslated_node-1];
      randomblength = randomblength - bl_array[untranslated_node-1];
      double reduce_length =   original_length - randomblength; 
      //cout<<"\nOriginal Length "<<original_length<<" randomblength "<<randomblength<<" reduce_length "<<reduce_length;
      
      
      
      /* change length of branch where inserted */
      t_weight[selected_edge][0] = randomblength;
            
      /* insert new internal node to attach new leaf */
      unsigned int i_n = t.new_node();
      t_weight[i_n][0] = reduce_length;
      //bl_array[currnodecount] = bl_array[currnodecount-1] + reduce_length;
  
      /* insert the new leaf */
      unsigned int l_n = t.new_node(); 
      //t_weight[l_n][0] = randomblength;/*NEED TO CHANGE THIS */
      //t_weight[l_n][0] = leafLength;
      t_name[l_n] = newleaf;
      parent_array[l_n] = i_n;
      
      /* Change LCA ARRAY if selected edge is the root i.e STEM CASE*/
      if(selected_edge == subtree_root_node && ADDoption[random_leaf] == STEM) {
        //cout<<"   CHANGING LCA !!";
        lca_array[subtree_root_node] = i_n;
      }
      
      /* Print single taxa family leaves */
      //if(single_taxa_family[random_leaf] == true) {
        //cout<<"\nSingle taxa family! selected edge = "<<selected_edge;
      //}
      /* remove the edge b/w selected node and parent */
      unsigned int current_parent = parent_array[selected_edge];//parent of selected edge
      //cout<<"\nCurrent parent "<<current_parent<<" insert branch "<<selected_edge;
      t.remove_edge(current_parent, selected_edge); 
            
      //cout<<"\nnot root selected";
      //cout<<"\nThe parent is : "<<parent;	      
     
     
          
      /* add an edge b/w new internal and selected node */
      t.add_edge(i_n, selected_edge);
      parent_array[selected_edge] = i_n;
          
      /* assign parent of new internal to old parent of selected */
      t.add_edge(i_n, current_parent);
      parent_array[i_n] = current_parent;
      
      
      double leafLength = 0;
      unsigned int parent = parent_array[i_n];
      //find branch-length of new leaf for ultra-metric tree
      /*Do a DFS from selected edge to leaf*/
      for (aw::Tree::iterator_dfs v=t.begin_dfs(i_n, parent),vEE=t.end_dfs(); v!=vEE; ++v) {
	unsigned int node = v.idx;
	//cout<<"\nNode "<<current_node<<" Branch Length: "<<t_weight[current_node][0];
	if (node == i_n)//don't add current nodes length 
	  continue;
	//cout<<"\n INSIDE DFS "<<node;
	leafLength += t_weight[node][0];
	if(t.is_leaf(node)){
	  //cout<<t_name[node];
	  break;	  
	}	  
      }        
      
       /* add an edge b/w new leaf and new internal */       
      t.add_edge(i_n, l_n);
      
      //cout<<"\nThe new calculated leaf length is "<<leafLength;
      
      //Weight of new leaf node
      t_weight[l_n][0] = leafLength; 
      
      /* Display the new number of edges and nodes */ 
      /*unsigned int edgec = t.edge_size();
      unsigned int nodec = t.node_size();
      cout<<"\nCurrent number of edges is "<<edgec;
      cout<<"\nCurrent number of nodes is "<<nodec;*/
      
      /* Output tree to verify */
      /*TREE_POSTORDER2(k,t){
	unsigned int currnode = k.idx;
	cout<<"\nWeight for node "<<currnode<<" "<<t_name[currnode];
	cout<<" is "<< t_weight[currnode][0];       
      } */  
      //int wait;
      //cin>>wait;
      delete [] bl_array; 
      ostringstream os;
      aw::tree2newick(os, t, t_name, t_weight);
      cout<<"\n"<<os.str()<<endl;
         
    }
    //unsigned int edgec = t.edge_size();
    //unsigned int  nodec = t.node_size();
    //cout<<"\nCurrent number of edges is "<<edgec;
    //cout<<"\nCurrent number of nodes is "<<nodec;
    
    
    double bl_temp =0;
    TREE_POSTORDER2(k,t) {
	unsigned int currnode = k.idx;    
	bl_temp += t_weight[currnode][0]; 
    }   
    cout<<"\nThe new total branch length is: "<<bl_temp;
    
    /*Write the tree to file */
    aw::tree2newick(ofs, t, t_name, t_weight);
    ofs<<endl;
   
    delete[] parent_array;
    delete[] leafindex; 
    delete[] added_leaf;
  }
 
  delete [] initial_parent; 
  delete [] translate_index; 
  
  cout<<"\n";
  return 0;
  
  
}  


int binarysearch(double *larray, int size, double key) {
  int first = 0;
  int last = size-1;
  if(size == 1) {
   return 0;
  }
  //cout<<"\nInside binary search!!";
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
