#if !defined(READINPUTFILE_INCLUDED)
#define READINPUTFILE_INCLUDED

void removeDuplicates(std::vector<int>& vec)
{
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

void removeDuplicates(std::vector<string>& vec)
{
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}


class readInputFile {

 private: 
  // temporary variables for this function
  int iS, iF, iT;
  string s,f,t;

  string w;
  double iW;
  //string Sso, Ssi, St;

  node *newnode;
  edge current;	    

  int newNodeIndexso;
  int newNodeIndexsi;
  int index;
    
  vector<node> ni;
  vector<node> selfInterNodesi;
  vector<edge> elisti;

  vector<node> temp;

  vector<string> old_s;
  vector<string> old_f;
  vector<string> old_we;  

  vector<int> old_iS;
  vector<int> old_iF;
  vector<double> old_iWE;  
  
  vector<string> UniqueStr;
  vector<int>    node_ids;

  vector<node> emptyn;
  vector<edge> emptye;
  vector<int>  emptyi;
  vector<string> emptys;

  Helper ihelper;

 public:
  readInputFile();
  ~readInputFile();

  bool areDigits   ( string _str );
  bool checkValues ( char *firstLine, int Ncols );

  void readFile    ( const char *file_network, int Ncols, bool Header, int skip, bool quite, int NalhpaNum );
 
  void enumerateNodeList( vector<node> nj, vector<edge> ej, vector<node> &Rnj, vector<edge> &Rej );
 
  vector<node> getNodeSet();
  vector<edge> getEdgeSet();

};

readInputFile::readInputFile(){}

/*
  Utility method returns true if a string contains only 
  integers from 0 to 9.
 */
bool readInputFile::areDigits( string _str ){

  for (int i=0; i<_str.length(); i++)
    {
      if( _str[i]>='0' && _str[i]<='9'){;}
      else
	return false;
          
    return true;
    }
}


/*
 Utility method to read the files head line, and the number 
 of data columns in the file.
*/
bool readInputFile::checkValues( char *firstLine, int Ncols){

  char *tokens = strtok(firstLine, "\t\r"); 
  vector<string> columns;     
  while(tokens != NULL){
    columns.push_back( string(tokens) );       
    tokens = strtok(NULL, "\t\r");     
  }

  if( areDigits(columns[0]) && areDigits(columns[1]) ){

    node_ids.push_back((int)atoi(columns[0].c_str()));
    node_ids.push_back((int)atoi(columns[1].c_str()));

    if( columns.size() == 2 ){
      old_iS.push_back((int)atoi(columns[0].c_str()));
      old_iF.push_back((int)atoi(columns[1].c_str()));
      old_iWE.push_back(1);
    } 

    if( columns.size() == 3 ){

      old_iS.push_back((int)atoi(columns[0].c_str()));
      old_iF.push_back((int)atoi(columns[1].c_str()));
      if( Ncols == 3 )
	old_iWE.push_back((double)atof(columns[2].c_str()));
      else
	old_iWE.push_back(1);

    }

    return false;

  } else {

    UniqueStr.push_back(columns[0]);
    UniqueStr.push_back(columns[1]);

    if( columns.size() == 2 ){
      old_s.push_back(columns[0]);
      old_f.push_back(columns[1]);
      old_we.push_back("1");
    }


    if( columns.size() == 3 ){
      old_s.push_back(columns[0]);
      old_f.push_back(columns[1]);
      if( Ncols == 3 )
	old_we.push_back(columns[2]);
      else
	old_we.push_back("1");
    }

    return true;

  }

  return true;

}


/*
 Read in network file which can be weighted or unweighted; 
 Ncols = 2 => unweighted
 Ncols = 3 => weighted 
 File input file can have a header line at the begining:
 Header = (True,False)
 Skip n number of lines starting at the begining. 
 */
void readInputFile::readFile(const char *file_network, int Ncols, bool Header, int skip, bool quite, int NalphaNum ){

  // Scan input file
  ifstream fscan(file_network, ios::in);
   
  if( fscan.is_open() != true ){
    cout << "> Error opening " << file_network << endl;
    exit(1);
  }

  if( quite ){ cout << "> Scanning network_file: " << file_network << endl; }

  // Print header 
  if(Header){
    string hh;
    getline(fscan,hh);
    if( quite ){ cout << "Header: " << hh << endl; }
  }

  // Skip lines from the begining
  if( skip > 0 ){
    string hh;
    for( int i=0; i<skip; i++ ){
      getline(fscan,hh);
    }
  }

  int count = 1;    
  char firstLine[256];
  bool stringValues=false;

  //read-in network data

  if( NalphaNum == -1 ){//test leading ids if alpha-numeric
    fscan.getline(firstLine, 256);
    stringValues = checkValues(firstLine,Ncols);
  } else {
   
    if( NalphaNum == 1 ){ stringValues = true;  }//dealing with string ids

    if( NalphaNum == 2 ){ stringValues = false; }//dealing with numeric ids
  
  }  

 
  if(stringValues){

    // Read in the datafile
    count = 1;
    if(Ncols == 2){
 
      while (fscan >> s >> f ) {	        // read friendship pair (s,f)
	if ( f.compare(s) < 0) { t = s; s = f; f = t; }	// guarantee s < f
	
	old_s.push_back(s);
	old_f.push_back(f);
	
	UniqueStr.push_back(s);
	UniqueStr.push_back(f);
	
	old_we.push_back("1");
      }
    } else {

      while (fscan >> s >> f >> w ) {	        // read friendship pair (s,f)
	if ( f.compare(s) < 0) { t = s; s = f; f = t; }	// guarantee s < f
      
	old_s.push_back(s);
	old_f.push_back(f);
      
	UniqueStr.push_back(s);
	UniqueStr.push_back(f);
	
	old_we.push_back(w);
      }	 
    }
     
    fscan.close();

    // Create list of unique node ids
    removeDuplicates(UniqueStr);

    // Remove duplicate edges
    for(int i=0; i < old_s.size(); ++i){
    
      for(int j=i+1; j < old_s.size(); ++j){
      
	//if( (old_s[j].compare(old_s[i]) == 0) && (old_f[j].compare(old_f[i]) == 0 ) )
	if( old_s[j] == old_s[i] && old_f[j] == old_f[i] )
	  { old_s[j] = "CDMSuite_NA"; old_f[j] = "CDMSuite_NA"; old_we[j] = "CDMSuite_NA"; }
      }
       
      if( (old_s[i].compare("CDMSuite_NA") != 0) && (old_f[i].compare("CDMSuite_NA") != 0) ){
	if(  atof(old_we[i].c_str()) == 0 ) old_we[i] = "1";
	
	int id_so = ihelper.rBinarySearch( UniqueStr, 0, UniqueStr.size(), old_s[i] );
	int id_si = ihelper.rBinarySearch( UniqueStr, 0, UniqueStr.size(), old_f[i] );

	elisti.push_back(edge( (id_so+1), (id_si+1), atof(old_we[i].c_str()), count++, (id_so+1), (id_si+1), old_s[i], old_f[i] ));
      } 
    }

    // Store node and edge lists
    count = 1;
    vector<int>  tally;
    for(int i=0; i<UniqueStr.size(); i++){
      tally.push_back( count );
      ni.push_back( node( count, 0.0, 0.0, (i+1), UniqueStr[i] ) );
      count++;
    }
     
    ni.insert(ni.begin(), node(-1, 0.0, 0.0, -1, "") );	      
 
    for(int i=0; i<ni.size(); i++){

      for(int j=0; j<elisti.size(); j++){
       
	if( elisti[j].IDso == ni[i].ID || elisti[j].IDsi == ni[i].ID )
	  ni[i].addEdge( elisti[j].so, elisti[j].si, elisti[j].we, elisti[j].key, elisti[j].IDso, elisti[j].IDsi, elisti[j].sIDso, elisti[j].sIDsi );
	
      }

    }

  } 



  if(!stringValues){

    // Read in the datafile
    count = 1;
    if(Ncols == 2){
 
      while (fscan >> iS >> iF ) {	        // read friendship pair (s,f)
	if ( iF < iS ) { iT = iS; iS = iF; iF = iT; }	// guarantee s < f
	
	old_iS.push_back(iS);
	old_iF.push_back(iF);
	
	old_iWE.push_back(1);

	node_ids.push_back( iS );
	node_ids.push_back( iF );

      }
    } else {

      while (fscan >> iS >> iF >> iW ) {	        // read friendship pair (s,f)
	if ( iF < iS ) { iT = iS; iS = iF; iF = iT; }	// guarantee s < f
	
	old_iS.push_back(iS);
	old_iF.push_back(iF);
      
	old_iWE.push_back(iW);

	node_ids.push_back( iS );
	node_ids.push_back( iF );

      }	 
    }
     
    fscan.close();

    
    // Create list of unique node ids
    removeDuplicates(node_ids);

    char tt1[256];
    char tt2[256];
    
    // Create list of unique node ids
    for(int i=0; i < old_iS.size(); ++i){
	
      for(int j=i+1; j < old_iS.size(); ++j){
	
	if( old_iS[j] == old_iS[i] && old_iF[j] == old_iF[i] )
	  {old_iS[j] = -1; old_iF[j] = -1; old_iWE[j] = -1;}
      }
      
      if(old_iS[i] != -1 && old_iF[i] != -1){
	if(old_iWE[i] == 0) old_iWE[i] = 1;
	
	int id_so = ihelper.rBinarySearch( node_ids, 0, node_ids.size(), old_iS[i] );
	int id_si = ihelper.rBinarySearch( node_ids, 0, node_ids.size(), old_iF[i] );

	sprintf(tt1,"%d",old_iS[i]);
	sprintf(tt2,"%d",old_iF[i]);
	
	elisti.push_back(edge( (id_so+1), (id_si+1), old_iWE[i], count++, (id_so+1), (id_si+1), tt1, tt2));

      } 
    }
    
    // Store node and edge lists
    count = 1;
    vector<int>  tally;
    for(int i=0; i<node_ids.size(); i++){
      tally.push_back( count );
      sprintf(tt1,"%d",node_ids[i]);
      ni.push_back( node( count, 0.0, 0.0, (i+1),tt1 ));
      count++;
    }
     
    ni.insert(ni.begin(), node(-1, 0.0, 0.0, -1, "") );	      
    
    for(int i=0; i<ni.size(); i++){
      
      for(int j=0; j<elisti.size(); j++){
       
	if( elisti[j].IDso == ni[i].ID || elisti[j].IDsi == ni[i].ID )
	  ni[i].addEdge( elisti[j].so, elisti[j].si, elisti[j].we, elisti[j].key, elisti[j].IDso, elisti[j].IDsi, elisti[j].sIDso, elisti[j].sIDsi );
       
      }

    }

  }

  //old_s.clear();
  //old_f.clear();
  //old_we.clear();
       
  return;

}

/*
 Renumerates the node and edge list after randomly removing nodes from the boot-strap method. 
 */
void readInputFile::enumerateNodeList(vector<node> nj, vector<edge> ej, vector<node> &Rnj, vector<edge> &Rej){
  
  ni.clear();
  elisti.clear();

  Rnj.clear();
  Rej.clear();

  node_ids.clear();

  int count = 1;
  vector<int> tally;
  tally.clear();
  for(int i=0; i<nj.size(); i++){  
    if( nj[i].c != -1 ){
      node_ids.push_back( nj[i].ID );
      tally.push_back( count );
      ni.push_back( node( count, 0.0, 0.0, nj[i].ID, nj[i].sID ) );          
      count++;
    }
  }

  count = 1;
  ni.insert(ni.begin(), node(-1, 0.0, 0.0, -1, "") );	      
  for(int i=0; i<ej.size(); i++){

    edge current = ej[i];
    
    if( current.so != -1 && current.si != -1 ){
      
      int index_so   = -1;
      int index_si   = -1;
      
      index_so       = ihelper.rBinarySearch( node_ids, 0, node_ids.size(), current.IDso );
      current.so     = tally[index_so];
            
      index_si       = ihelper.rBinarySearch( node_ids, 0, node_ids.size(), current.IDsi );
      current.si     = tally[index_si];

      if( index_so != -1 && index_si != -1 )
	elisti.push_back( edge(current.so, current.si, current.Globalwe, count++, current.IDso, current.IDsi, current.sIDso, current.sIDsi) );

    }
  }

  for(int i=0; i<ni.size(); i++){
    
    for(int j=0; j<elisti.size(); j++){
      
    if( elisti[j].IDso == ni[i].ID || elisti[j].IDsi == ni[i].ID )
      ni[i].addEdge( elisti[j].so, elisti[j].si, elisti[j].Globalwe, elisti[j].key, elisti[j].IDso, elisti[j].IDsi, elisti[j].sIDso, elisti[j].sIDsi );
      
    }

  }

  Rnj = ni;
  Rej = elisti;
  
  ni.clear();
  elisti.clear();
  


}

readInputFile::~readInputFile(){}

//return the new node list after reading the input file
vector<node> readInputFile::getNodeSet(){ return ni; }

//return the new edge list after reading the input file
vector<edge> readInputFile::getEdgeSet(){ return elisti; }


#endif
