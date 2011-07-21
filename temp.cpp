
  /* CREATE REPLICATES OF NEW  TREE */
  for(int k=0; k<replicates; k++) {
    int currnodecount = initial_currnodecount;
    //cout<<"Initial curr node count is: "<<currnodecount;
    double total_blength = total_blength_initial;
    cout<<"\n\nReplicate number : "<<k+1;
    cout<<"\nThe total branch length is "<<total_blength;
    double* bl_array = new double[total_nodes];
    int* parent_array = new int[total_nodes];
    unsigned int* translate_index = new unsigned int[total_nodes];
    for(int i=0; i<currnodecount; i++) {
      translate_index[i] = translate_index_initial[i];
      parent_array[i] = parent_array_initial[i];
      bl_array[i] = bl_array_initial[i];
    }
  
    /* Declare tree, labels & weights */
    aw::Tree t = initial_t;
    aw::idx2name t_name = initial_t_name;
    aw::idx2weight_double t_weight = initial_t_weight;
    
    //aw::Tree t;
    //aw::idx2name t_name;
    //aw::idx2weight_double t_weight;
    //aw::stream2tree(ifs, t, t_name, t_weight);
    
    /* Initialise all leaves as not added and copy leaf indices */
    for(int j=0; j<leafcount; j++) {
      leafindex[j] = j;
      added_leaf[j] = NOTADDED;
    }
  
    //int edgec = t.edge_size();
    //int nodec = t.node_size();
    //cout<<"\nBEfore addition Number of edges is "<<edgec;
    //cout<<"\nBefore addition Number of nodes is "<<nodec;
    int templc = leafcount;
    //cout<<"\nNumber of leaves to be added "<<templc;
  
    for(int j=0; j<leafcount; j++) {
    
    
      int edgec = t.edge_size();
      int nodec = t.node_size();
      //cout<<"\nInside Number of edges is "<<edgec;
      //cout<<" Inside Number of nodes is "<<nodec;
      /* Pick a random leaf */
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
      //cout<<"Leaf number to be added "<<newleaf;
    
      //int temp = leafindex[templc-1];
      //leafindex[templc-1] = leafindex[rnum] ;
      //leafindex[rnum] = temp;
      //templc--; 
      //for(int j=0; j<leafcount; j++) {
      //    cout<<"\nleaf index "<<leafindex[j];
      //    cout<<"  added "<<added_leaf[j];
      //}
    
    
    
      /* Pick a random length to select branch */
      //cout<<"\nTotal blength "<<total_blength;
      double randomblength;
      
      //cout<<"\nRandomblength "<<randomblength;  
      //double randomblength = fmod((double)rand() , (double)total_blength );
    
      /* Find the node to insert leaf based on branch lengths */
      unsigned int untranslated_node;
      do {
        randomblength = (double)rand() * (double)total_blength / (double)RAND_MAX; 
        //cout<<"\nTemp total branch length"<<total_blength;
        //cout<<" temp random branch length"<<randomblength;
        untranslated_node = binarysearch( bl_array, currnodecount-1, randomblength);  
      } while(translate_index[untranslated_node] == t.root);
      
      double original_length;
      //if(untranslated_node!=0)
      //original_length = bl_array[untranslated_node] - bl_array[untranslated_node-1];
      //else 
    
      /* obtain the individual lengths by subtracting the cumulative lengths */
      original_length = bl_array[untranslated_node] - bl_array[untranslated_node-1];
      randomblength = randomblength - bl_array[untranslated_node-1];
    
      double reduce_length =   original_length - randomblength; 
      //cout<<" reduce length"<<reduce_length;
      unsigned int insert_branch = translate_index[untranslated_node]; 
      //cout<<"\nnormalized Randomblength "<<randomblength;
      //cout<<"\nBranch where to insert is: "<<insert_branch; 
      //cout<<"\nReduce Length "<<reduce_length<<" Randomblength "<<randomblength<<" original length "<<original_length;
      //cout<<"\nnew curr node count "<<currnodecount;
    
      /* change length of branch where inserted */
      t_weight[insert_branch][0] = randomblength;
   
      /* insert new internal node to attach new leaf */
      unsigned int i_n = t.new_node();
      t_weight[i_n][0] = reduce_length;    
      //bl_array[currnodecount] = bl_array[currnodecount-1] + reduce_length;
  
      /* insert the new leaf */
      unsigned int l_n = t.new_node(); 
      t_weight[l_n][0] = randomblength; 
    
      //bl_array[currnodecount] = bl_array[currnodecount-1] + randomblength;
      //translate_index[currnodecount] = l_n;
      //currnodecount++;
      //currleafcount++;
    
     
      /* remove the edge b/w selected node and parent */
      unsigned int current_parent = parent_array[insert_branch];
      //cout<<"\nCurrent parent "<<current_parent<<" insert branch "<<insert_branch;
      t.remove_edge(current_parent, insert_branch); 
      
      
      //cout<<"\nnot root selected";
      //cout<<"\nThe parent is : "<<parent;	      
     
      /* add an edge b/w new leaf and new internal */ 
      t_name[l_n] = newleaf;
      t.add_edge(i_n, l_n);
      parent_array[l_n] = i_n;
    
      /* add an edge b/w new internal and selected node */
      t.add_edge(i_n, insert_branch);
      parent_array[insert_branch] = i_n;
      t.add_edge(i_n, current_parent);
    
      /* assign parent of new internal to old parent of selected */
      parent_array[i_n] = current_parent; 
      edgec = t.edge_size();
      nodec = t.node_size();
      
      total_blength = 0;
      currnodecount = 0;
    
      TREE_POSTORDER2(k,t){
	unsigned int currnode = k.idx;
	//cout<<"\nWeight for node "<<currnode<<" "<<t_name[currnode];
	//cout<<" is "<< t_weight[currnode][0];    
	total_blength += t_weight[currnode][0]; 
	translate_index[currnodecount] = currnode;
	if(currnode!=t.root)
	  bl_array[currnodecount++] = total_blength;//this array stores cumulative branch lengths & hence is sorted.
	//if(t.is_leaf(k.idx)){
        //leafn++;
	//}
	parent_array[currnode] = k.parent;     
      }   
      //currleafcount = leafn;
    
    
      /* print the branch lengths */
      
      /*for(int i=0; i<currnodecount-1; i++){
	int translate = translate_index[i]; 
	cout<<"\nbranch number "<<translate;
	if(i!=0)cout<<" length "<<bl_array[i] - bl_array[i-1];
	else cout<<" length "<<bl_array[i];
	cout<<" cumul length "<<bl_array[i];
	}
      
    
	//cout<<"\n";
	//TREE_POSTORDER2(k,t) {
	//cout<<k.idx<<" ";
	//} 
	 Print out tree */
      //ostringstream os;
      //aw::tree2newick(os, t, t_name, t_weight);
      //cout<<"\nOutput : "<<os.str()<<endl;
    
      //int wait;
      //cin>>wait;
    } 
    //edgec = t.edge_size();
    //nodec = t.node_size();
    //cout<<"\nAfter addition Number of edges is "<<edgec;
    //cout<<"\nAfter addition Number of nodes is "<<nodec;
    delete [] parent_array;
    delete [] bl_array;
    aw::tree2newick(ofs, t, t_name, t_weight);
    ofs<<endl;
    cout<<"\nThe new total branch length is: "<<total_blength;
  }
