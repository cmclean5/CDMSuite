#if !defined(ENRICHMENTSTUDIES_INCLUDED)
#define ENRICHMENTSTUDIES_INCLUDED

class enrichmentStudies {
  
 private:
  vector<vector<double> > pv_values;

  vector<string> study_IDs;
  vector<string> study_names;
  vector<int>    entrez_ids; 
  vector<int>    comsize;
  vector<int>    comsize_rnd;
  vector<int>    studies_tally;

  vector<vector<int> > studies;

  vector<vector<double> > permute;  
  //vector<vector<double> > permute_upper;
  vector<vector<double> > permute_lower;

  vector<vector<double> > permute_ll;  
  vector<vector<double> > permute_ll2;  

  vector<vector<double> > permute_log;  
  vector<vector<double> > permute_log2;  

  vector<double> permute_anno_log;  
  vector<double> permute_anno_log2;

  vector<int> comids;
  vector<int> coms;
  vector<int> coms_rnd;
  vector<string> ids1;
  vector<string> ids2;
  fstream *fileout20;
  fstream *fileout21;
  fstream *fileout22;

  typedef std::pair<double,int> mypair;

  void CalculateFDR1( vector<mypair> order_pv, int M, double level, double &fdr, double &pv, int &cutoff, vector<double> &_vals );
  void CalculateFDR2( vector<mypair> order_pv, int M, double level, double &fdr, double &pv, int &cutoff, vector<double> &_vals );
  double  _min(double,double);
  int    _min(int,int);
  //double binomial_co(int,int);
  //double factorial( double );
  double prob_overlap( double N, double na, double nb, double nab );
  void print_message(const char*);

  struct sort_pred{
    bool operator()(const std::pair<double,int> &l, const std::pair<double,int> &r){
      return l.first < r.first;
    }
  };


 public:
  enrichmentStudies();
  ~enrichmentStudies();
  void calculateDEPC(const char*, const char*, const char*);
  void overlapAnnotationAandB( int _N, const char* File1, const char* File2, const char* ext );
  void overlapAnnotationAandBinNetwork( const char* File1, const char* File2, const char* ext );
  void overlapAnnotationAandCommunities( const char* File1, const char* ext );
  void calculate_Disease_Anno_Co( const char*File1, const char* File2, const char* ext );
  void hypergeometricTest_rndComs();
  void hypergeometricTest_rndAnnoList();
  void hypergeometricTest2();
  void hypergeometricTest();
  void hypergeometricTest( fstream* fout1, fstream* fout2);
  void permutation        ( double );
  void permutation_coms   ( double );
  void permutation_comsize( double );
};

enrichmentStudies::enrichmentStudies(){}

enrichmentStudies::~enrichmentStudies(){}


//--- FDR1 -> Benjamini and Hochberg FDR (BH)
void enrichmentStudies::CalculateFDR1( vector<mypair> order_pv, int M, double level, double &fdr, double &pv, int &cutoff, vector<double> &_vals ){

  if( M <= 0 ) return;  
  int i = M;
  
  fdr    = level * (double)1/M;
  pv     = order_pv[0].first;
  cutoff = 1;

  for( vector<mypair>::iterator it=order_pv.end()-1; it != order_pv.begin()-1; it--, i--){  

    mypair xx       = (*it);
    double _pv      = xx.first;
    double pv_test  = level * (double) i/ (double) M;

    _vals.push_back( pv_test );

    if( pv_test < _pv ){
      fdr    = pv_test;
      pv     = _pv;
      cutoff = i;
      //return;
    } 
  }
      

}

//--- FDR2 -> Benjamini and Liu FDR (BL)
void enrichmentStudies::CalculateFDR2( vector<mypair> order_pv, int M, double level, double &fdr, double &pv, int &cutoff, vector<double> &_vals ){

  
  if( M <= 0 ) return;
  int i = 1;
  
  fdr    = level;
  pv     = order_pv[M-1].first;
  cutoff = M;

  for( vector<mypair>::iterator it=order_pv.begin(); it != order_pv.end(); it++, i++){  
    
    mypair xx       = (*it);
    double _pv      = xx.first;
    double pv_test  = level * M / ( (M + 1 - i) * (M + 1 - i) );

    if( level < pv_test ) pv_test = level;

    _vals.push_back( pv_test );

    if( pv_test > _pv ){
      fdr    = pv_test;
      pv     = _pv;
      cutoff = i;
      //return;
    }
    
    
  }
      

}

double enrichmentStudies::_min( double a, double b ){ return a<=b ? a : b; }

int    enrichmentStudies::_min( int a, int b ){ return a<=b ? a : b; }

/*
double enrichmentStudies::factorial( double n ){

  long double fact;
  //double fact;
  fact=1;

  while (n>0){    
    fact=fact*n;    
    n=n-1;    
  }    

  return fact;

}
*/

/*
double enrichmentStudies::binomial_co(int n, int k) {
    
    if (k==0 || n==k) return 1;
    if (n<=0 || k<0 || n<k) return 0;
  
    int k1=_min(k,n-k);
    int k2=n-k1;
    double fact=k2+1;
    for (double i=k1;i>1.;--i)
      fact *= (k2+i)/i;
    return fact;
  
}
*/

void enrichmentStudies::print_message(const char *message){cout << string(message) << endl;}


//calculateDEPC (Disease Enrichment Per Cluster)
//File1 = community file, Flie2 = DiseaseIDs.csv, File3 = EntrezIDs.csv
void enrichmentStudies::calculateDEPC(const char *File1, const char *File2, const char *ext){

  string ids1i;
  int ids2i, ids3i;

  int comi1, comi2;
  int com_max = 1;

  char comments[256];

  ifstream filein;
  filein.open(File1);
  filein.getline(comments, 256);
  print_message(comments);


  char line[256];
  //--- Read community file
  while( filein.getline(line,256) ){

    const char *del = " = \r";
    char *tokens    = strtok(line, del);
    vector<string> columns;
    while(tokens !=NULL){
      columns.push_back( string(tokens) ); 
      tokens = strtok(NULL, del);
    }
    if( columns.size() < 2 ) continue;

    comi1 = atoi(columns[0].c_str()) ;
    comi2 = atoi(columns[1].c_str()) ;

    comids.push_back  ( comi1 );
    coms.push_back    ( comi2 );
    coms_rnd.push_back( -1 );

    if( comi2 > com_max )
      com_max = comi2;

  }
  

  filein.close();  
  
  cout << "com_max " << com_max << endl;
  comsize.resize(com_max);
  comsize_rnd.resize(com_max);
  int com   = 1;
  int _size = 0;
  int tally = 0;

  while(com<(com_max+1)){
    _size = 0;
    for(int i=0; i<comids.size(); i++){
      if(coms[i] == com ){ _size++; }
    }
    if(_size != 0)
      comsize[tally++] = _size;
    com++;
  }


  //--- combinations
  filein.open(File2);
  filein.getline(comments, 256);
  print_message(comments);


  //--- Read flatfile
  while( filein.getline(line,256) ){

    const char *del = " \t \r";
    char *tokens    = strtok(line, del);
    vector<string> columns;
    while(tokens !=NULL){
      columns.push_back( string(tokens) ); 
      tokens = strtok(NULL, del);
    }
    if( columns.size() < 2 ) continue;

    study_IDs.push_back  ( columns[0].c_str() ) ;
    ids1.push_back       ( columns[0].c_str() ) ;
    ids2.push_back       ( columns[1].c_str() ) ;
    entrez_ids.push_back ( atoi(columns[2].c_str()) ) ;

  }
  
  filein.close();  

  //---unique list of disease IDs & names
  std::sort(study_IDs.begin(),study_IDs.end());
  study_IDs.erase(std::unique(study_IDs.begin(),study_IDs.end()),study_IDs.end());

  for(int s=0; s<study_IDs.size(); s++){  
    for( int i=0; i<ids1.size(); i++){
      if( ids1[i].compare(study_IDs[s]) == 0 ){
	bool found = false;	
	for(int n=0; n<study_names.size(); n++){
	  if( ids2[i].compare(study_names[n]) == 0 ) found = true;
	}
	if(!found)
	  study_names.push_back( ids2[i] );
      }
    }
  }

  int no_ids = 0;

  studies_tally.resize(study_names.size());
  studies.resize(comids.size());
  permute.resize(comids.size());
  //permute_upper.resize(comids.size());
  permute_lower.resize(comids.size());

  for (int i = 0; i < comids.size(); i++){
    studies[i].resize(study_IDs.size());
    permute[i].resize(study_IDs.size());
    //permute_upper[i].resize(study_IDs.size());
    permute_lower[i].resize(study_IDs.size());
  }

  for(int k=0; k<entrez_ids.size(); k++){

    int search = entrez_ids[k];

    int ind1 = -1;
    int ind2 = -1;
 
    for( int i=0; i< comids.size(); i++){
      if( search == comids[i] ){
	ind1 = k;
	ind2 = i;
      }
    }  

    if( ind1 != -1 ){

      for(int s=0; s<study_IDs.size(); s++){    

	if( ids1[ind1].compare(study_IDs[s]) == 0 )
	  studies[ind2][s] = 1;
      }


    } else {;
    }

  }


  //---table entrez_id verses ext(Disease,function,compartment)
  for(int s=0; s<studies.size(); s++){
    for(int i=0; i<studies[0].size(); i++){
      studies_tally[i] += studies[s][i];
    }
  }
  
  int tot = 0;
  for(int s=0; s<study_IDs.size(); s++){ 
    tot += studies_tally[s];
  }

  cout << "" << endl;
  cout << "Tol: " << tot << endl;

  //---size of network
  double N = comids.size();

  //p-values (if network a single community)
  char buffer  [250];
  char buffer2 [250];
 
  //--- p-value calculations
  sprintf(buffer ,"OUT/p_values_%s.csv",ext);
  sprintf(buffer2,"OUT/enrichment_%s.csv",ext);
  fstream *fileout5 = new fstream(buffer,ios_base::out);
  fstream *fileout6 = new fstream(buffer2,ios_base::out);

  cout << "Network: "        << (int)comids.size() << " Communities: " << com_max << endl;

  cout         << "numbers for entire network" << endl;
  (*fileout5)  << "Community" << "\t" << "Size" << "\t";
  (*fileout6)  << "Community" << "\t" << "Size" << "\t";
  for(int p=0; p<studies_tally.size(); p++){
    cout         << study_names[p] << " " << studies_tally[p] << " "; 
    (*fileout5)  << study_names[p] << "\t";
    (*fileout6)  << study_names[p] << "\t";
  }

  cout << "" << endl;
  cout << "--------" << endl;

  pv_values.resize(comsize.size());
  for (int i = 0; i < comsize.size(); i++)
    pv_values[i].resize(study_IDs.size());

  //---calculate p-values per annotation per cluster.
  hypergeometricTest( fileout5, fileout6 );
 
  //---permutation of annotations in list to test p-value robustness
  //--- p-value calculations
  sprintf(buffer ,"OUT/permute_p_values_%s.csv",ext);
  fstream *fileout15 = new fstream(buffer,ios_base::out);

  if(true){

    double P=1000;
    double J=2;

    double Bon1 = 0.05/(double)(comsize.size()*studies_tally.size());
    double Bon2 = 0.01/(double)(comsize.size()*studies_tally.size());

    char buffer3  [250];
    char buffer4  [250];

    for( int p=0; p<P; p++){
      
      permutation( (p+1) );//permute ids in the annotation list
    
      hypergeometricTest_rndAnnoList();

      //permutation_coms( (p+1) );//permute ids pairs between communities
      //hypergeometricTest_rndComs();

    }

    (*fileout15)  << "Community" << "\t" << "Size" << "\t";

    for(int p=0; p<studies_tally.size(); p++){
      //sprintf(buffer3,"Bonferroni 0.05 significance level (%e)",Bon1);
      //sprintf(buffer4,"Bonferroni 0.01 significance level (%e)",Bon2);
      //(*fileout15)  << study_names[p] << "\t" << "p-value" << "\t" << "Permutated p_value <= (%)" <<"\t" << buffer3 << "\t" << buffer4 << "\t";
      (*fileout15)  << study_names[p] << "\t" << "p-value" << "\t" << "Permutated p_value <= (%)" <<"\t" << "Bonferroni 0.05 significance level" << "\t" << "Bonferroni 0.01 significance level" << "\t";
    }
    
   
    //--- loop over all communities
    for( int m=0; m<comsize.size(); m++){

      (*fileout15) << "" << endl;
      (*fileout15) << (m+1) << "\t" << comsize[m] << "\t";

      //--- loop over each Disease type
      //The significance level sqrt( p-value * (1-p-value)/P )
      for( int f=0; f<studies_tally.size(); f++ ){
	string str1 = "-";
	string str2 = "-";
	if( pv_values[m][f] <= Bon1 ){ str1 = "**";  }
	if( pv_values[m][f] <= Bon2 ){ str2 = "***"; }
	(*fileout15) << "\t" << pv_values[m][f] << "\t" << permute_lower[m][f]/P * 100 << "\t" << str1 << "\t" << str2 << "\t";
      }  
    }

  }
    
  fileout15->close();


  
  //--- Calculate FDR
  vector<vector<double> > fdr1_values;
  vector<vector<double> > fdr2_values;
  vector<vector<double> > PV_values;
  vector<vector<int> >    com_values;
  fdr1_values.resize(comsize.size());
  fdr2_values.resize(comsize.size());
  PV_values.resize(comsize.size());
  com_values.resize(comsize.size());
  for (int i = 0; i < comsize.size(); i++){
    fdr1_values[i].resize(study_IDs.size());
    fdr2_values[i].resize(study_IDs.size());
    PV_values[i].resize(study_IDs.size());
    com_values[i].resize(study_IDs.size());
  }

  sprintf(buffer,"OUT/p_sig_tests_%s.csv",ext);
  fstream *fileout11 = new fstream(buffer,ios_base::out);
  for(int f=0; f<studies_tally.size(); f++){    
    (*fileout11) << study_names[f] << "\t" << "" << "\t" << "" << "\t" << "" << "\t" << "" << "\t" << "" << "\t" << "" << "\t" << "";      
  }
  (*fileout11) << "" << endl;
  for(int f=0; f<studies_tally.size(); f++){
    (*fileout11) << "Community" << "\t" << "Size" << "\t" << "p-values" << "\t" << "{p-values}" << "\t" << "Rank" << "\t" << "FDR (BH)" << "\t" << "FDR (BL)" << "\t";
  }
  (*fileout11) << "" << endl;
  
  sprintf(buffer,"OUT/p_sig_tests_summary_%s.csv",ext);
  fstream *fileout12 = new fstream(buffer,ios_base::out);
  (*fileout12) << "DISEASE" << "\t" << "p-value" << "\t" << "{p-value}" << "\t" << "FDR (BH)" << "\t" << "LEVEL" << "\t" << "p-value" << "\t" << "{p-value}" << "\t" << "FDR (BL)" << "\t" << "LEVEL" << endl; 


  for(int f=0; f<studies_tally.size(); f++){
    
    vector<mypair> pv_sorted;  
    (*fileout12) << study_names[f] << "\t";
    for(int m=0; m<comsize.size();m++){
      pv_sorted.push_back( mypair(pv_values[m][f],m) );
    }
  
    std::sort( pv_sorted.begin(), pv_sorted.end(), sort_pred() );
 
    double FDR = 0.0; double PV = 0.0; int LEVEL = 0; double STAT_LEVEL = 0.05; vector<double> vals; vector<double> vals2;
    CalculateFDR1( pv_sorted, comsize.size(), STAT_LEVEL, FDR, PV, LEVEL, vals );
    (*fileout12) << PV << "\t" << PV*(double)(comsize.size()/LEVEL) << "\t" << FDR << "\t" << LEVEL << "\t";
    
    FDR = 0.0; PV = 0.0; LEVEL = 0;
    CalculateFDR2( pv_sorted, comsize.size(), STAT_LEVEL, FDR, PV, LEVEL, vals2 );
    (*fileout12) << PV << "\t" << PV*(double)(comsize.size()/LEVEL) << "\t" << FDR << "\t" << LEVEL << endl;

    for(int m=0; m<comsize.size();m++){
      fdr1_values[m][f] = vals[m];
      fdr2_values[m][f] = vals2[m];
      PV_values[m][f]   = pv_sorted[m].first;
      com_values[m][f]  = pv_sorted[m].second;
    }

  }

  int i=1;
  int M=comsize.size()-1;
  for(int m=0; m<comsize.size();m++,M--){
    for(int f=0; f<studies_tally.size(); f++){
      (*fileout11) << (com_values[m][f]+1) << "\t" << comsize[com_values[m][f]] << "\t" << PV_values[m][f] << "\t" <<  PV_values[m][f] * (double)( (comsize.size())/i ) << "\t" << i << "\t" << fdr1_values[M][f] << "\t" << fdr2_values[m][f] << "\t";
    }
    i++;
    (*fileout11) << "" << endl;
  }

  

  fileout11->close();
  fileout12->close();
  

  //cout << "DONE" << endl;

}

double enrichmentStudies::prob_overlap( double N, double na, double nb, double nab ){

  double result = 0.0;

  double num = factln( (int)na ) + factln( (int)(N-na) ) + factln( (int)nb ) + factln( (int)(N-nb) );
  
  double dem = factln( (int)N )  + factln( (int)(na-nab) ) + factln( (int)nab ) + factln( (int)(N - na - nb + nab) ) + factln( (int)(nb - nab) );

  result = exp( num - dem ) ;


  return result;
  
}


//N, File1 = data/flatfile_MASC_family.csv, File2 = data/flatfile_MASC_disease.csv 
void enrichmentStudies::overlapAnnotationAandB( int _N, const char* File1, const char* File2, const char* ext ){

  //needed to find overlap with annotation A
  vector<string> n1;
  vector<string> n1_ID;
  vector<int>    n1_ids;
  vector<int>    n1_ids_i;
  vector<string> n1_name;
  vector<string> n1_type;

  //needed to find overlap with annotation B
  vector<string> n2;
  vector<string> n2_ID;
  vector<int>    n2_ids;
  vector<int>    n2_ids_j;
  vector<string> n2_name;
  vector<string> n2_type;

  char comments[256];
  
  ifstream filein;
  filein.open(File1);
  filein.getline(comments, 256);
  //print_message(comments);

  char line[256];
  //--- Read community file
  while( filein.getline(line,256) ){

    const char *del = "\t\r";
    char *tokens    = strtok(line, del);
    vector<string> columns;
    while(tokens !=NULL){
      columns.push_back( string(tokens) ); 
      tokens = strtok(NULL, del);
    }
    if( columns.size() < 2 ) continue;

    n1.push_back      ( columns[0].c_str()  );
    n1_ID.push_back   ( columns[0].c_str()  );
    n1_name.push_back ( columns[1].c_str()  );
    n1_ids.push_back  ( atoi(columns[2].c_str()) );
    n1_ids_i.push_back( atoi(columns[2].c_str()) );

  }

  filein.close();  

  //---unique list of n1 names 
  std::sort(n1_ID.begin(),n1_ID.end());
  n1_ID.erase(std::unique(n1_ID.begin(),n1_ID.end()),n1_ID.end());

  for(int s=0; s<n1_ID.size(); s++){  
    for( int i=0; i<n1_ids.size(); i++){
      if( n1[i].compare(n1_ID[s]) == 0 ){
	bool found = false;	
	for(int n=0; n<n1_type.size(); n++){
	  if( n1_name[i].compare(n1_type[n]) == 0 ) found = true;
	}
	if(!found)
	  n1_type.push_back( n1_name[i] );
      }
    }
  }

  filein.open(File2);
  filein.getline(comments, 256);


  //--- Read community file
  while( filein.getline(line,256) ){

    const char *del = "\t\r";
    char *tokens    = strtok(line, del);
    vector<string> columns;
    while(tokens !=NULL){
      columns.push_back( string(tokens) ); 
      tokens = strtok(NULL, del);
    }
    if( columns.size() < 2 ) continue;

    n2.push_back      ( columns[0].c_str() );
    n2_ID.push_back   ( columns[0].c_str() );
    n2_name.push_back ( columns[1].c_str() );
    n2_ids.push_back  ( atoi(columns[2].c_str()) );
    n2_ids_j.push_back( atoi(columns[2].c_str()) );    
  }

  filein.close();  

  
  //---unique list of n2 names 
  std::sort(n2_ID.begin(),n2_ID.end());
  n2_ID.erase(std::unique(n2_ID.begin(),n2_ID.end()),n2_ID.end());

  for(int s=0; s<n2_ID.size(); s++){  
    for( int i=0; i<n2_ids.size(); i++){
      if( n2[i].compare(n2_ID[s]) == 0 ){
	bool found = false;	
	for(int n=0; n<n2_type.size(); n++){
	  if( n2_name[i].compare(n2_type[n]) == 0 ) found = true;
	}
	if(!found)
	  n2_type.push_back( n2_name[i] );
      }
    }
  }

  
  //---unique list of n1 IDs 
  std::sort(n1_ids_i.begin(),n1_ids_i.end());
  n1_ids_i.erase(std::unique(n1_ids_i.begin(),n1_ids_i.end()),n1_ids_i.end());

  //---unique list of n2 IDs 
  std::sort(n2_ids_j.begin(),n2_ids_j.end());
  n2_ids_j.erase(std::unique(n2_ids_j.begin(),n2_ids_j.end()),n2_ids_j.end());

  vector<int> all_ids;
  for( int i=0; i<n1_ids_i.size(); i++)
    all_ids.push_back( n1_ids_i[i] );

  for( int i=0; i<n2_ids_j.size(); i++)
    all_ids.push_back( n2_ids_j[i] );

  //---unique list of all IDs 
  std::sort(all_ids.begin(),all_ids.end());
  all_ids.erase(std::unique(all_ids.begin(),all_ids.end()),all_ids.end());


  //overlap matrix: rows (n1), columns (n2)
  vector<vector<double> > overlap;
  overlap.resize(n1_type.size());

  vector<vector<double> > muab;
  muab.resize(n1_type.size());

  for(int i=0; i<n1_type.size(); i++){
    overlap[i].resize(n2_type.size());
    muab[i].resize(n2_type.size());
  }

  vector<double> na;
  na.resize(n1_type.size());

  vector<double> nb;
  nb.resize(n2_type.size());

  int ind_c [na.size()]; 
  int ind_r [nb.size()]; 

  //---print out overlap info
  char buffer  [250];
  sprintf(buffer ,"OUT/overlap_%s.csv",ext);
  fstream *fileout = new fstream(buffer,ios_base::out);
      
  
  double N = _N;

  for(int _na=0; _na<na.size(); _na++ ) ind_c[_na] = 0;
  for(int _nb=0; _nb<nb.size(); _nb++ ) ind_r[_nb] = 0;
	
  for(int x=0; x<n1_type.size(); x++){
    for(int y=0; y<n2_type.size(); y++){
      overlap[x][y] = 0.0;
    }
  }

  //---loop over all ids in network and find overlap between annotation n1 and annotation type n2
  for( int i=0; i<all_ids.size(); i++ ){

    for(int _na=0; _na<na.size(); _na++ ) ind_c[_na] = 0;
    for(int _nb=0; _nb<nb.size(); _nb++ ) ind_r[_nb] = 0;
	
    //loop over all annotation type 'n1'
    for( int _na=0; _na<n1_ids.size(); _na++){
	  
      if( (n1_ids[_na] == all_ids[i]) ){
	

	for(int s=0; s<n1_type.size(); s++){
	  int temp_ind_c = -1;
	  
	  if( n1_type[s].compare(n1_name[_na]) == 0 ){ ind_c[s] = 1; temp_ind_c = s;}
	  
	  if( temp_ind_c != -1 ) { na[temp_ind_c]++; }
	  
	}

      }
    }

    //loop over all annotations: annotation type 'n2'
    for( int _nb=0; _nb<n2_ids.size(); _nb++){
      
      if( (n2_ids[_nb] == all_ids[i]) ){
	
	for(int s=0; s<n2_type.size(); s++){
	  int temp_ind_r = -1;
	      
	  if( n2_type[s].compare(n2_name[_nb]) == 0 ){ ind_r[s] = 1; temp_ind_r = s; }
	  
	  if( temp_ind_r != -1 ) { nb[temp_ind_r]++; }
	}
	    
      }
    }
    

    //entrez ID in network and overlap between disease and annotation
    for(int x=0; x<na.size(); x++){
      for(int y=0; y<nb.size(); y++){
	if( ind_c[x] == 1 && ind_r[y] == 1 )
	  overlap[x][y] += 1;
      }
    }	
    
  }
  
 
  (*fileout) << "N: " << N << " unique ids in File 1: " << n1_ids_i.size() << " unique ids in File 2: " << n2_ids_j.size() << " all ids: " << all_ids.size() << endl;
  (*fileout) << "annotation 1" << "\t" << "n1" << "\t" << "annotation 2" << "\t" << "n2" << "\t" << "actual overlap" << "\t" << "p-value" << endl;
 
  //---reset and calculate muab
  for(int _na=0; _na<n1_type.size(); _na++){
    for(int _nb=0; _nb<n2_type.size(); _nb++){
      muab[_na][_nb] = 0.0;
      muab[_na][_nb] = prob_overlap( (double)N, (double)na[_na], (double)nb[_nb], (double)overlap[_na][_nb] );  
    }
  }	

  //---calculate prob
  for(int _na=0; _na<n1_type.size(); _na++){
    for(int _nb=0; _nb<n2_type.size(); _nb++){
      double prob_tot = 0.0;
      for(int k=0; k<=overlap[_na][_nb]; k++){
	double prob = prob_overlap( (double)N, (double)na[_na], (double)nb[_nb], (double)k );
	if( prob <= muab[_na][_nb] ) prob_tot += prob;
      }
	  
      //---rounding error
      if( prob_tot > 1 ) prob_tot = 1.0;
      
     
      (*fileout) << n1_type[_na] << "\t" << na[_na] << "\t" << n2_type[_nb] << "\t" << nb[_nb] << "\t" << 
	overlap[_na][_nb] << "\t" << prob_tot << endl;
  
    }
  }

  
  fileout->close();
  
 
}

//File1 = data/flatfile_MASC_family.csv, File2 = data/flatfile_MASC_disease.csv 
void enrichmentStudies::overlapAnnotationAandBinNetwork( const char* File1, const char* File2, const char* ext ){

  //needed to find overlap with annotation A
  vector<string> n1;
  vector<string> n1_ID;
  vector<int>    n1_ids;
  vector<string> n1_name;
  vector<string> n1_type;

  //needed to find overlap with annotation B
  vector<string> n2;
  vector<string> n2_ID;
  vector<int>    n2_ids;
  vector<string> n2_name;
  vector<string> n2_type;

  char comments[256];
  
  ifstream filein;
  filein.open(File1);
  filein.getline(comments, 256);
  //print_message(comments);

  char line[256];
  //--- Read community file
  while( filein.getline(line,256) ){

    const char *del = "\t\r";
    char *tokens    = strtok(line, del);
    vector<string> columns;
    while(tokens !=NULL){
      columns.push_back( string(tokens) ); 
      tokens = strtok(NULL, del);
    }
    if( columns.size() < 2 ) continue;

    n1.push_back      ( columns[0].c_str()  );
    n1_ID.push_back   ( columns[0].c_str()  );
    n1_name.push_back ( columns[1].c_str()  );
    n1_ids.push_back  ( atoi(columns[2].c_str()) );

  }

  filein.close();  

  //---unique list of n1 names 
  std::sort(n1_ID.begin(),n1_ID.end());
  n1_ID.erase(std::unique(n1_ID.begin(),n1_ID.end()),n1_ID.end());

  for(int s=0; s<n1_ID.size(); s++){  
    for( int i=0; i<n1_ids.size(); i++){
      if( n1[i].compare(n1_ID[s]) == 0 ){
	bool found = false;	
	for(int n=0; n<n1_type.size(); n++){
	  if( n1_name[i].compare(n1_type[n]) == 0 ) found = true;
	}
	if(!found)
	  n1_type.push_back( n1_name[i] );
      }
    }
  }

  filein.open(File2);
  filein.getline(comments, 256);


  //--- Read community file
  while( filein.getline(line,256) ){

    const char *del = "\t\r";
    char *tokens    = strtok(line, del);
    vector<string> columns;
    while(tokens !=NULL){
      columns.push_back( string(tokens) ); 
      tokens = strtok(NULL, del);
    }
    if( columns.size() < 2 ) continue;

    n2.push_back      ( columns[0].c_str() );
    n2_ID.push_back   ( columns[0].c_str() );
    n2_name.push_back ( columns[1].c_str() );
    n2_ids.push_back  ( atoi(columns[2].c_str()) );
    
  }

  filein.close();  

  //---unique list of n2 names 
  std::sort(n2_ID.begin(),n2_ID.end());
  n2_ID.erase(std::unique(n2_ID.begin(),n2_ID.end()),n2_ID.end());

  for(int s=0; s<n2_ID.size(); s++){  
    for( int i=0; i<n2_ids.size(); i++){
      if( n2[i].compare(n2_ID[s]) == 0 ){
	bool found = false;	
	for(int n=0; n<n2_type.size(); n++){
	  if( n2_name[i].compare(n2_type[n]) == 0 ) found = true;
	}
	if(!found)
	  n2_type.push_back( n2_name[i] );
      }
    }
  }

  
  //overlap matrix: rows (n1), columns (n2)
  vector<vector<double> > overlap;
  overlap.resize(n1_type.size());

  vector<vector<double> > muab;
  muab.resize(n1_type.size());

  for(int i=0; i<n1_type.size(); i++){
    overlap[i].resize(n2_type.size());
    muab[i].resize(n2_type.size());
  }

  vector<double> na;
  na.resize(n1_type.size());

  vector<double> nb;
  nb.resize(n2_type.size());

  int ind_c [na.size()]; 
  int ind_r [nb.size()]; 

  //---print out overlap info
  char buffer  [250];
  sprintf(buffer ,"OUT/overlap_%s.csv",ext);
  fstream *fileout = new fstream(buffer,ios_base::out);
      
  double N = comids.size();

  for(int _na=0; _na<na.size(); _na++ ) ind_c[_na] = 0;
  for(int _nb=0; _nb<nb.size(); _nb++ ) ind_r[_nb] = 0;
	
  for(int x=0; x<n1_type.size(); x++){
    for(int y=0; y<n2_type.size(); y++){
      overlap[x][y] = 0.0;
    }
  }

  //---loop over all ids in network and find overlap between annotation n1 and annotation type n2
  for( int i=0; i<comids.size(); i++ ){

    for(int _na=0; _na<na.size(); _na++ ) ind_c[_na] = 0;
    for(int _nb=0; _nb<nb.size(); _nb++ ) ind_r[_nb] = 0;
	
    //loop over all annotation type 'n1'
    for( int _na=0; _na<n1_ids.size(); _na++){
	  
      if( (n1_ids[_na] == comids[i]) ){
	

	for(int s=0; s<n1_type.size(); s++){
	  int temp_ind_c = -1;
	  
	  if( n1_type[s].compare(n1_name[_na]) == 0 ){ ind_c[s] = 1; temp_ind_c = s;}
	  
	  if( temp_ind_c != -1 ) { na[temp_ind_c]++; }
	  
	}

      }
    }

    //loop over all annotations: annotation type 'n2'
    for( int _nb=0; _nb<n2_ids.size(); _nb++){
      
      if( (n2_ids[_nb] == comids[i]) ){
	
	for(int s=0; s<n2_type.size(); s++){
	  int temp_ind_r = -1;
	      
	  if( n2_type[s].compare(n2_name[_nb]) == 0 ){ ind_r[s] = 1; temp_ind_r = s; }
	  
	  if( temp_ind_r != -1 ) { nb[temp_ind_r]++; }
	}
	    
      }
    }
    

    //entrez ID in network and overlap between disease and annotation
    for(int x=0; x<na.size(); x++){
      for(int y=0; y<nb.size(); y++){
	if( ind_c[x] == 1 && ind_r[y] == 1 )
	  overlap[x][y] += 1;
      }
    }	
    
  }

   
  (*fileout) << "N: " << N << endl;
  (*fileout) << "annotation 1" << "\t" << "n1" << "\t" << "annotation 2" << "\t" << "n2" << "\t" << "actual overlap" << "\t" << "p-value" << endl;
 
  //---reset and calculate muab
  for(int _na=0; _na<n1_type.size(); _na++){
    for(int _nb=0; _nb<n2_type.size(); _nb++){
      muab[_na][_nb] = 0.0;
      muab[_na][_nb] = prob_overlap( (double)N, (double)na[_na], (double)nb[_nb], (double)overlap[_na][_nb] );  
    }
  }	

  //---calculate prob
  for(int _na=0; _na<n1_type.size(); _na++){
    for(int _nb=0; _nb<n2_type.size(); _nb++){
      double prob_tot = 0.0;
      for(int k=0; k<=overlap[_na][_nb]; k++){
	double prob = prob_overlap( (double)N, (double)na[_na], (double)nb[_nb], (double)k );
	if( prob <= muab[_na][_nb] ) prob_tot += prob;
      }
	  
      //---rounding error
      if( prob_tot > 1 ) prob_tot = 1.0;
      
      //cout << n1_type[_na] << "\t" << na[_na] << "\t" << n2_type[_nb] << "\t" << nb[_nb] << "\t" << 
      //overlap[_na][_nb] << "\t" << prob_tot << endl;

      (*fileout) << n1_type[_na] << "\t" << na[_na] << "\t" << n2_type[_nb] << "\t" << nb[_nb] << "\t" << 
	overlap[_na][_nb] << "\t" << prob_tot << endl;
  
    }
  }

  
  fileout->close();
  
  //cout << "DONE" << endl;
  
}

//File1 = data/flatfile_MASC_family.csv, File2 = data/flatfile_MASC_disease.csv 
void enrichmentStudies::overlapAnnotationAandCommunities( const char* File1, const char* ext ){

  //needed to find overlap with annotation A
  vector<string> n1;
  vector<string> n1_ID;
  vector<int>    n1_ids;
  vector<string> n1_name;
  vector<string> n1_type;

  //needed to find overlap with annotation B
  vector<string> n2;
  vector<string> n2_ID;
  vector<int>    n2_ids;
  vector<int>    n2_ids_j;
  vector<string> n2_name;
  vector<string> n2_type;

  char line[256];

  //--- annotation type 1 is the clustering  
  int com_max = 0;
  for( int i=0; i<comids.size(); i++ ){
    n1_ids.push_back( comids[i] );
    sprintf( line, "community_%d",coms[i] );
    n1_name.push_back( line );
    if( coms[i] > com_max )
      com_max = coms[i];
  }

  for( int i=0; i<com_max; i++){
    sprintf(line, "community_%d", (i+1) );
    n1_type.push_back( line ); 
  }


  char comments[256];
  
  ifstream filein;
  filein.open(File1);
  filein.getline(comments, 256);
  //print_message(comments);

  //--- Read community file
  while( filein.getline(line,256) ){

    const char *del = "\t\r";
    char *tokens    = strtok(line, del);
    vector<string> columns;
    while(tokens !=NULL){
      columns.push_back( string(tokens) ); 
      tokens = strtok(NULL, del);
    }
    if( columns.size() < 2 ) continue;

    n2.push_back      ( columns[0].c_str()  );
    n2_ID.push_back   ( columns[0].c_str()  );
    n2_name.push_back ( columns[1].c_str()  );
    n2_ids.push_back  ( atoi(columns[2].c_str()) );
    n2_ids_j.push_back( atoi(columns[2].c_str()) ); 

  }

  filein.close();  

  //---unique list of n1 names 
  std::sort(n2_ID.begin(),n2_ID.end());
  n2_ID.erase(std::unique(n2_ID.begin(),n2_ID.end()),n2_ID.end());

  for(int s=0; s<n2_ID.size(); s++){  
    for( int i=0; i<n2_ids.size(); i++){
      if( n2[i].compare(n2_ID[s]) == 0 ){
	bool found = false;	
	for(int n=0; n<n2_type.size(); n++){
	  if( n2_name[i].compare(n2_type[n]) == 0 ) found = true;
	}
	if(!found)
	  n2_type.push_back( n2_name[i] );
      }
    }
  }

  //overlap matrix: rows (n1), columns (n2)
  vector<vector<double> > overlap;
  overlap.resize(n1_type.size());

  vector<vector<double> > muab;
  muab.resize(n1_type.size());

  for(int i=0; i<n1_type.size(); i++){
    overlap[i].resize(n2_type.size());
    muab[i].resize(n2_type.size());
  }

  vector<double> na;
  na.resize(n1_type.size());

  vector<double> nb;
  nb.resize(n2_type.size());

  int ind_c [na.size()]; 
  int ind_r [nb.size()]; 

  //---print out overlap info
  char buffer  [250];
  sprintf(buffer ,"OUT/overlap_%s.csv",ext);
  fstream *fileout = new fstream(buffer,ios_base::out);
      
  
  double N = comids.size();

  for(int _na=0; _na<na.size(); _na++ ) ind_c[_na] = 0;
  for(int _nb=0; _nb<nb.size(); _nb++ ) ind_r[_nb] = 0;
	
  for(int x=0; x<n1_type.size(); x++){
    for(int y=0; y<n2_type.size(); y++){
      overlap[x][y] = 0.0;
    }
  }

  //---loop over all ids in network and find overlap between annotation n1 and annotation type n2
  for( int i=0; i<comids.size(); i++ ){

    for(int _na=0; _na<na.size(); _na++ ) ind_c[_na] = 0;
    for(int _nb=0; _nb<nb.size(); _nb++ ) ind_r[_nb] = 0;
	
    //loop over all annotation type 'n1'
    for( int _na=0; _na<n1_ids.size(); _na++){
	  
      if( n1_ids[_na] == comids[i] ){
	

	for(int s=0; s<n1_type.size(); s++){
	  int temp_ind_c = -1;
	  
	  if( n1_type[s].compare(n1_name[_na]) == 0 ){ ind_c[s] = 1; temp_ind_c = s;}
	  
	  if( temp_ind_c != -1 ) { na[temp_ind_c]++; }
	  
	}

      }
    }

    //loop over all annotations: annotation type 'n2'
    for( int _nb=0; _nb<n2_ids.size(); _nb++){
      
      if( n2_ids[_nb] == comids[i] ){
	
	for(int s=0; s<n2_type.size(); s++){
	  int temp_ind_r = -1;
	      
	  if( n2_type[s].compare(n2_name[_nb]) == 0 ){ ind_r[s] = 1; temp_ind_r = s; }
	  
	  if( temp_ind_r != -1 ) { nb[temp_ind_r]++; }
	}
	    
      }
    }
    

    //entrez ID in network and overlap between disease and annotation
    for(int x=0; x<na.size(); x++){
      for(int y=0; y<nb.size(); y++){
	if( ind_c[x] == 1 && ind_r[y] == 1 )
	  overlap[x][y] += 1;
      }
    }	
    
  }
  
  
  (*fileout) << "N: " << N << endl;
  (*fileout) << "annotation 1" << "\t" << "n1" << "\t" << "annotation 2" << "\t" << "n2" << "\t" << "actual overlap" << "\t" << "p-value" << endl;
 
  //---reset and calculate muab
  for(int _na=0; _na<n1_type.size(); _na++){
    for(int _nb=0; _nb<n2_type.size(); _nb++){
      muab[_na][_nb] = 0.0;
      muab[_na][_nb] = prob_overlap( (double)N, (double)na[_na], (double)nb[_nb], (double)overlap[_na][_nb] );  
    }
  }	

  //---calculate prob
  for(int _na=0; _na<n1_type.size(); _na++){
    for(int _nb=0; _nb<n2_type.size(); _nb++){
      double prob_tot = 0.0;
      for(int k=0; k<=overlap[_na][_nb]; k++){
	double prob = prob_overlap( (double)N, (double)na[_na], (double)nb[_nb], (double)k );
	if( prob <= muab[_na][_nb] ) prob_tot += prob;
      }
	  
      //---rounding error
      if( prob_tot > 1 ) prob_tot = 1.0;
      
     
      (*fileout) << n1_type[_na] << "\t" << na[_na] << "\t" << n2_type[_nb] << "\t" << nb[_nb] << "\t" << 
	overlap[_na][_nb] << "\t" << prob_tot << endl;
  
    }
  }

  
  fileout->close();
  
 
}


//File1 = flatfile, File2 = p_value_network  
void enrichmentStudies::calculate_Disease_Anno_Co( const char* File1, const char* File2, const char* ext ){
  
   //needed to find overlap with annotation A
  vector<string> n1;
  vector<string> n1_ID;
  vector<int>    n1_ids;
  vector<string> n1_name;
  vector<string> n1_type;

  //needed to find overlap with annotation B
  vector<string> n2;
  vector<string> n2_ID;
  vector<int>    n2_ids;
  vector<string> n2_name;
  vector<string> n2_type;

  char comments[256];
  
  ifstream filein;
  filein.open(File1);
  filein.getline(comments, 256);
  //print_message(comments);

  char line[256];
  //--- Read community file
  while( filein.getline(line,256) ){

    const char *del = "\t\r";
    char *tokens    = strtok(line, del);
    vector<string> columns;
    while(tokens !=NULL){
      columns.push_back( string(tokens) ); 
      tokens = strtok(NULL, del);
    }
    if( columns.size() < 2 ) continue;

    n1.push_back      ( columns[0].c_str()  );
    n1_ID.push_back   ( columns[0].c_str()  );
    n1_name.push_back ( columns[1].c_str()  );
    n1_ids.push_back  ( atoi(columns[2].c_str()) );

  }

  filein.close();  

  //---unique list of n1 names 
  std::sort(n1_ID.begin(),n1_ID.end());
  n1_ID.erase(std::unique(n1_ID.begin(),n1_ID.end()),n1_ID.end());

  for(int s=0; s<n1_ID.size(); s++){  
    for( int i=0; i<n1_ids.size(); i++){
      if( n1[i].compare(n1_ID[s]) == 0 ){
	bool found = false;	
	for(int n=0; n<n1_type.size(); n++){
	  if( n1_name[i].compare(n1_type[n]) == 0 ) found = true;
	}
	if(!found)
	  n1_type.push_back( n1_name[i] );
      }
    }
  }

  filein.open(File2);
  filein.getline(comments, 256);


  //--- Read community file
  while( filein.getline(line,256) ){

    const char *del = "\t\r";
    char *tokens    = strtok(line, del);
    vector<string> columns;
    while(tokens !=NULL){
      columns.push_back( string(tokens) ); 
      tokens = strtok(NULL, del);
    }
    if( columns.size() < 2 ) continue;

    n2.push_back      ( columns[0].c_str() );
    n2_ID.push_back   ( columns[0].c_str() );
    n2_name.push_back ( columns[1].c_str() );
    n2_ids.push_back  ( atoi(columns[2].c_str()) );
    
  }

  filein.close();  

  
  //---unique list of n2 names 
  std::sort(n2_ID.begin(),n2_ID.end());
  n2_ID.erase(std::unique(n2_ID.begin(),n2_ID.end()),n2_ID.end());

  for(int s=0; s<n2_ID.size(); s++){  
    for( int i=0; i<n2_ids.size(); i++){
      if( n2[i].compare(n2_ID[s]) == 0 ){
	bool found = false;	
	for(int n=0; n<n2_type.size(); n++){
	  if( n2_name[i].compare(n2_type[n]) == 0 ) found = true;
	}
	if(!found)
	  n2_type.push_back( n2_name[i] );
      }
    }
  }

  
  //overlap matrix: rows (n1), columns (n2)
  vector<vector<double> > overlap;
  overlap.resize(n1_type.size());

  vector<vector<double> > muab;
  muab.resize(n1_type.size());

  for(int i=0; i<n1_type.size(); i++){
    overlap[i].resize(n2_type.size());
    muab[i].resize(n2_type.size());
  }

  vector<double> na;
  na.resize(n1_type.size());

  vector<double> nb;
  nb.resize(n2_type.size());

  int ind_c [na.size()]; 
  int ind_r [nb.size()]; 

  //---print out overlap info
  char buffer  [250];
  sprintf(buffer ,"OUT/overlap_%s.csv",ext);
  fstream *fileout = new fstream(buffer,ios_base::out);
      
  
  //double N = comids.size();
      
  //---loop over each community
  for(int n=0; n<comsize.size(); n++){

    if( comsize[n] > 2 ){

      double N   = (double)comsize[n];

      for(int _na=0; _na<na.size(); _na++ ) na[_na] = 0;
      for(int _nb=0; _nb<nb.size(); _nb++ ) nb[_nb] = 0;

      for(int x=0; x<n1_type.size(); x++){
	for(int y=0; y<n2_type.size(); y++){
	  overlap[x][y] = 0.0;
	}
      }

      //---loop over all ids in community (comids) and find overlap between disease & annotation, or annotation type 'a' and annotation type 'b'
      for( int i=0; i<comids.size(); i++ ){

	for(int _na=0; _na<na.size(); _na++ ) ind_c[_na] = 0;
	for(int _nb=0; _nb<nb.size(); _nb++ ) ind_r[_nb] = 0;
	
	//loop over all disease: annotation type 'a'
	for( int _na=0; _na<n1_ids.size(); _na++){
	  
	  if( (n1_ids[_na] == comids[i]) && ((n+1) == coms[i]) ){

	    for(int s=0; s<n1_type.size(); s++){
	      int temp_ind_c = -1;
	      
	      if( n1_type[s].compare(n1_name[_na]) == 0 ){ ind_c[s] = 1; temp_ind_c = s;}
	    
	      if( temp_ind_c != -1 ) { na[temp_ind_c]++; }

	    }

	  }
	}

	//loop over all annotations: annotation type 'b'
	for( int _nb=0; _nb<n2_ids.size(); _nb++){
	    
	  if( (n2_ids[_nb] == comids[i]) && ((n+1) == coms[i]) ){	 

	    for(int s=0; s<n2_type.size(); s++){
	      int temp_ind_r = -1;
		
	      if( n2_type[s].compare(n2_name[_nb]) == 0 ){ ind_r[s] = 1; temp_ind_r = s; }
	      
	      if( temp_ind_r != -1 ) { nb[temp_ind_r]++; }
	    }
	    
	  }
	}


	//entrez ID in network and overlap between disease and annotation
	for(int _na=0; _na<na.size(); _na++){
	  for(int _nb=0; _nb<nb.size(); _nb++){
	    if( ind_c[_na] == 1 && ind_r[_nb] == 1 && (n+1) == coms[i] )	    
	      overlap[_na][_nb] += 1;
	  }
	}	

      }
        
      (*fileout) << "community" << "\t" << (n+1) << "\t" << "size" << "\t" << comsize[n] << endl;
      (*fileout) << "" << "\t";      
      for(int p=0; p<n2_type.size(); p++){
	(*fileout)  << n2_type[p] << "\t";
      }
 
      (*fileout) << "" << endl;

      //---reset and calculate muab
      for(int _na=0; _na<n1_type.size(); _na++){
	for(int _nb=0; _nb<n2_type.size(); _nb++){
	  muab[_na][_nb] = 0.0;
	  muab[_na][_nb] = prob_overlap( (double)N, (double)na[_na], (double)nb[_nb], (double)overlap[_na][_nb] );  
	}
      }	

      //---calculate prob
      for(int _na=0; _na<n1_type.size(); _na++){
	(*fileout) << n1_type[_na] << "\t";
	for(int _nb=0; _nb<n2_type.size(); _nb++){
	  double prob_tot = 0.0;
	  for(int k=0; k<=overlap[_na][_nb];k++){
	    double prob = prob_overlap( (double)N, (double)na[_na], (double)nb[_nb], (double)k );
	    if( prob <= muab[_na][_nb] ) prob_tot += prob;
	  }
	  
	  //---rounding error
	  if( prob_tot > 1 ) prob_tot = 1.0;

	  (*fileout) << prob_tot << "\t";
	    
	}
	(*fileout) << "" << endl;
      }

      (*fileout) << "" << endl;

    }
  }
 
  fileout->close();
  
 
}


void enrichmentStudies::hypergeometricTest( fstream* fout1, fstream* fout2){

   //---size of network
  double N = comids.size();

  //--- loop over all communities
  for( int m=0; m<comsize.size(); m++){

    (*fout1) << "" << endl;
    (*fout1) << (m+1) << "\t" << comsize[m] << "\t";

    (*fout2) << "" << endl;
    (*fout2) << (m+1) << "\t" << comsize[m] << "\t";

    

    //double n_choose_m = binomial_co( N, comsize[m] );
    double n_choose_m = bico( (int)N, (int)comsize[m] );
  
    //--- loop over each Disease type
    for( int f=0; f<studies_tally.size(); f++ ){

      double _f       = studies_tally[f];  
  
      double p_value  = 0.0;
    
      double en       = 0.0;

      int tally       = 0;

      //--- loop over all genes in the mth cluster which shares the fth Disease type
      for(int i=0; i<comids.size(); i++){

	if( coms[i] == (m+1) ){

	  if( studies[i][f] == 1 ){
      
	    //double nmf_choose_mmi = binomial_co( (N-_f), (comsize[m]-tally) );
	    //double f_choose_i     = binomial_co( _f, tally );
	    double nmf_choose_mmi = bico( (int)(N-_f), (int)(comsize[m]-tally) );
	    double f_choose_i     = bico( (int)_f, (int)tally );

	    if(isfinite( (f_choose_i * nmf_choose_mmi)/n_choose_m) == 1 )
	      p_value = p_value + (f_choose_i * nmf_choose_mmi)/n_choose_m;

	    tally++;
	    en      = tally;
	
	  }

	}
      }
    
      permute[m][f]         = 0.0;
      //permute_upper[m][f]   = 0.0;    
      permute_lower[m][f]   = 0.0;    


      pv_values[m][f] = (double)(1.0 - p_value);
      (*fout1) << (1.0-p_value)      << "\t";
      (*fout2) << en/(double)comsize[m] << "\t";      
    
    }
      
  }

  fout1->close();
  fout2->close();

}

void enrichmentStudies::hypergeometricTest_rndAnnoList(){

  double norm     = (double)comsize.size()*studies_tally.size();

  //---size of network
  double N = comids.size();

  double temp_p_values[studies_tally.size()];
  for(int i=0; i<studies_tally.size(); i++)
    temp_p_values[i] = 0.0;


  //--- loop over all communities
  for( int m=0; m<comsize.size(); m++){

    double n_choose_m = bico( (int)N, (int)comsize[m] );
      
    //--- loop over each Disease type
    for( int f=0; f<studies_tally.size(); f++ ){

      double _f       = studies_tally[f];  
  
      double p_value  = 0.0;
    
      int tally       = 0;

      //--- loop over all genes in the mth cluster which shares the fth Disease type
      for(int i=0; i<comids.size(); i++){

	if( coms[i] == (m+1) ){
	
	  if( studies[i][f] == 1 ){
      
	    double nmf_choose_mmi = bico( (int)(N-_f), (int)(comsize[m]-tally) );
	    double f_choose_i     = bico( (int)_f, (int)tally );

	    if(isfinite( (f_choose_i * nmf_choose_mmi)/n_choose_m) == 1 )
	      p_value = p_value + (f_choose_i * nmf_choose_mmi)/n_choose_m;

	    tally++;
	
	  }

	}
      }
    
      //correct for multiple-testing
      //Test every permute[m][f] value against every pv_values[m][f] value
      for( int M=0; M<comsize.size(); M++){
	for( int F=0; F<studies_tally.size(); F++ ){
	  //if( (double)(1.0 - p_value) >= pv_values[M][F] ) 
	  //{ permute_upper[M][F]   = permute_upper[M][F] + 1/norm; }
	  if( (double)(1.0 - p_value) <= pv_values[M][F] )
	    { permute_lower[M][F]   = permute_lower[M][F] + 1/norm; }
	}
      }
      
      permute[m][f]         = permute[m][f] + (double)(1.0 - p_value);
    
    }
      
  }

}

void enrichmentStudies::hypergeometricTest_rndComs(){

  double norm     = (double)comsize.size()*studies_tally.size();

  //---size of network
  double N = comids.size();

  double temp_p_values[studies_tally.size()];
  for(int i=0; i<studies_tally.size(); i++)
    temp_p_values[i] = 0.0;

  //--- loop over all communities
  for( int m=0; m<comsize.size(); m++){

    double n_choose_m = bico( (int)N, (int)comsize_rnd[m] );
  
    //--- loop over each Disease type
    for( int f=0; f<studies_tally.size(); f++ ){

      double _f       = studies_tally[f];  
  
      double p_value  = 0.0;
    
      int tally       = 0;

      //--- loop over all genes in the mth cluster which shares the fth Disease type
      for(int i=0; i<comids.size(); i++){

	if ( coms_rnd[i] == (m+1) ){

	  if( studies[i][f] == 1 ){
      
	    double nmf_choose_mmi = bico( (int)(N-_f), (int)(comsize_rnd[m]-tally) );
	    double f_choose_i     = bico( (int)_f, (int)tally );

	    if(isfinite( (f_choose_i * nmf_choose_mmi)/n_choose_m) == 1 )
	      p_value = p_value + (f_choose_i * nmf_choose_mmi)/n_choose_m;

	    tally++;
	
	  }

	}
      }
    
      //correct for multiple-testing
      //Test every permute[m][f] value against every pv_values[m][f] value
      for( int M=0; M<comsize.size(); M++){
	for( int F=0; F<studies_tally.size(); F++ ){
	  //if( (double)(1.0 - p_value) >= pv_values[M][F] )
	  //{ permute_upper[M][F]   = permute_upper[M][F] + 1/norm; }
	  if( (double)(1.0 - p_value) <= pv_values[M][F] )
	    { permute_lower[M][F]   = permute_lower[M][F] + 1/norm; }
	}
      }


      permute[m][f]         = permute[m][f] + (double)(1.0 - p_value);

    }
      
  }

}

void enrichmentStudies::hypergeometricTest2(){

   
  for(int p=0; p<20;p++){

    permutation_coms((double)p);

    //---size of network
    double N = comids.size();

    //--- loop over all communities
    for( int m=0; m<comsize.size(); m++){
      
      //double n_choose_m = binomial_co( N, comsize[m] );
      double n_choose_m = bico( (int)N, (int)comsize[m] );
  
      //--- loop over each Disease type
      for( int f=0; f<studies_tally.size(); f++ ){

	double _f       = studies_tally[f];  
	
	double p_value  = 0.0;
	
	int tally       = 0;

	//--- loop over all genes in the mth cluster which shares the fth Disease type
	for(int i=0; i<comids.size(); i++){

	  if( coms_rnd[i] == (m+1) ){

	    if( studies[i][f] == 1 ){
      
	      //double nmf_choose_mmi = binomial_co( (N-_f), (comsize[m]-tally) );
	      //double f_choose_i     = binomial_co( _f, tally );
	      double nmf_choose_mmi = bico( (int)(N-_f), (int)(comsize[m]-tally) );
	      double f_choose_i     = bico( (int)_f, (int)tally );
	      
	      if(isfinite( (f_choose_i * nmf_choose_mmi)/n_choose_m) == 1 )
		p_value = p_value + (f_choose_i * nmf_choose_mmi)/n_choose_m;
	      
	      tally++;
	      
	    }
	  }
	}

 
      }      
    }   
  }   


}



//---randomise the distribution of gene id's for each annotation
void enrichmentStudies::permutation( double tt ){

  time_t timer;
  struct tm y2k;
  double seconds;

  y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
  y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;

  time(&timer);  /* get current time; same as: timer = time(NULL)  */

  seconds = difftime(timer,mktime(&y2k)) * tt;


  //---create a new random number generator in current scope
  Ran _rand2(0);
  //--- Initialize random seed to the system clock. 
  _rand2.setSeed(seconds);

 
  //---reset the permutation matrix 
  for(int f=0; f<studies_tally.size(); f++){
    for(int i=0; i<comids.size(); i++){
      studies[i][f]      = 0;
    }
  }
  

  for(int f=0; f<studies_tally.size(); f++){  

    double F = studies_tally[f];

    int    N = comids.size()-1;

    int sample_size = 0;

    while( sample_size < F ){

      int rnd_ind = (int)floor( _rand2.doub() * (N) );

      if( studies[rnd_ind][f] == 0 ){
	studies[rnd_ind][f] = 1;
	sample_size++;
      }

    }

  }

 
}

//---randomise the distribution of gene id's for each annotation
void enrichmentStudies::permutation_comsize( double tt ){

  time_t timer;
  struct tm y2k;
  double seconds;

  y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
  y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;

  time(&timer);  /* get current time; same as: timer = time(NULL)  */

  seconds = difftime(timer,mktime(&y2k)) * tt;


  //---create a new random number generator in current scope
  Ran _rand2(0);
  //--- Initialize random seed to the system clock. 
  _rand2.setSeed(seconds);

  int C = comsize.size();

  int sample_size = 0;
 
  //--reset comsize
  for(int m=0; m<comsize.size(); m++){
    comsize_rnd[m] = 0;
  }

  //---reset the permutation matrix 
  for(int i=0; i<comids.size(); i++){
    coms_rnd[i] = -1;  
   }

  for(int i=0; i<comids.size(); i++){
    
    int rnd_C  = -1;
    
    while ( rnd_C <= 0 || rnd_C > C ){
      rnd_C  = (int)floor( _rand2.doub() * (C) ); 
    }

    coms_rnd[i] = rnd_C;
    comsize_rnd[rnd_C-1]++;

  }


}

//---randomise the distribution of gene id's for each annotation
void enrichmentStudies::permutation_coms( double tt ){

  time_t timer;
  struct tm y2k;
  double seconds;

  y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
  y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;

  time(&timer);  /* get current time; same as: timer = time(NULL)  */

  seconds = difftime(timer,mktime(&y2k)) * tt;


  //---create a new random number generator in current scope
  Ran _rand2(0);
  //--- Initialize random seed to the system clock. 
  _rand2.setSeed(seconds);

  Ran _rand3(0);
  _rand3.setSeed(seconds+10);

  //int C = comsize.size();

  int N = comids.size()-1;

  //vector<int> sample_size;
  //sample_size.resize( comsize.size() );
 
  //--set comsize
  for(int m=0; m<comsize.size(); m++){
    comsize_rnd[m] = comsize[m];
  //  sample_size[m] = comsize[m];
  }

  int sample = 2000;//N;

  //---reset the permutation matrix 
  for(int i=0; i<comids.size(); i++)
    coms_rnd[i] = coms[i];      
  
  while( sample > 0 ){

          int rnd_ind1 = (int)floor( _rand2.doub() * (N) );
          int rnd_ind2 = (int)floor( _rand3.doub() * (N) );

	  int temp_c   = coms_rnd[rnd_ind1];
	  coms_rnd[rnd_ind1] = coms_rnd[rnd_ind2];
	  coms_rnd[rnd_ind2] = temp_c;
	  //coms_rnd[rnd_ind2] = coms[rnd_ind1];
	  
	  sample--;
  }

}



/*
  for(int i=0; i<comids.size(); i++){
    
    bool get_next  = false;
    int tries      = 0;

    while(!get_next && tries < 1000 ){

      int rnd_C  = -1;

      while ( rnd_C <= 0 || rnd_C > C ){
	rnd_C  = (int)floor( _rand2.doub() * (C) ); 
      }

      if( sample_size[rnd_C-1] > 0 ){
	coms_rnd[i] = rnd_C;
	sample_size[rnd_C-1]--;	
	get_next = true;
      }

      tries++;
    }

  }

}
*/

#endif
