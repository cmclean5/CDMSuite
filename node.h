#if !defined(NODE_INCLUDED)
#define NODE_INCLUDED

/*
 The node class which contains: a unique key; it's 
 inputted ID from the file; a vector edges connected to this 
 node, and it's corresponding degree. 
 */
class node {

 public:
  double             w; //weight
  double             d; //distance
  int                k; //key
  int               ID; //Input ID from file 
  string           sID; //Input ID from file 
  int                v; //visited
  int                c; //community
  int           degree; //degree
  bool selfInteraction; //contains self-interacting edges

  string          name;

  node();
  ~node();
  node(int, double, double);
  node(int, double, double, int);
  node(int, double, double, int, string);
  node(int, double, double, string);

  void visited(int);
  void addEdge(edge newEdge);
  void addEdge(int s, int f, double w, int _k, int I, int F);
  void addEdge(int s, int f, double w, int _k, int I, int F, string sI, string sF);
  void addEdge(int s, int f, double w, int _k, string sI, string sF);
  void removeEdge(int edgeKey);
  void updateNodeCommunity(int newCom);
  vector<edge> getEdges();
  int getDegree();
  void setDegree( int newDegree );
  void print();
  void print(fstream* fout );
  void printEdges();
  node Clone();

 private:
  vector<edge> edges;//list of edges connected to this node.  

};

node::node(){k=-1; w=-1; d=-1; v=0; c=0; ID=0; sID=""; degree=0; selfInteraction = false; edges.clear(); name=""; }

node::node(int key, double we, double dist){ k = key; w = we; d = dist; v=0; c=0; ID=0; sID=""; degree=0; selfInteraction = false; edges.clear(); name=""; }

node::node(int key, double we, double dist, int _id){ k = key; w = we; d = dist; v=0; c=0; ID=_id; sID=""; degree=0; selfInteraction = false; edges.clear(); name=""; }

node::node(int key, double we, double dist, int _id, string _sid){ k = key; w = we; d = dist; v=0; c=0; ID=_id; sID=_sid; degree=0; selfInteraction = false; edges.clear(); name=""; }

node::node(int key, double we, double dist, string _sid){ k = key; w = we; d = dist; v=0; c=0; ID=0; sID=_sid; degree=0; selfInteraction = false; edges.clear(); name=""; }

node::~node(){}

void node::visited(int vis){ v = vis; }

void node::addEdge(edge newEdge){ 
  
  int eso    = newEdge.so;
  int esi    = newEdge.si;
  double ew  = newEdge.we;
  int ek     = newEdge.key;

  edge enew = edge(eso, esi, ew, ek);

  if(eso != esi){
    edges.push_back(enew);
  } else {
    if(!selfInteraction){
      edges.push_back(enew);
      selfInteraction = true;
    }
  }
}

void node::addEdge(int s, int f, double _w, int _k, int I, int F){ 
  
  edge enew = edge(s, f, _w, _k, I, F);

  if(s != f){
    edges.push_back(enew); 
  } else {
    if(!selfInteraction){
      edges.push_back(enew); 
      selfInteraction = true;
    }
  }

}

void node::addEdge(int s, int f, double _w, int _k, int I, int F, string sI, string sF){ 
  
  edge enew = edge(s, f, _w, _k, I, F, sI, sF);

  if(s != f){
    edges.push_back(enew); 
  } else {
    if(!selfInteraction){
      edges.push_back(enew); 
      selfInteraction = true;
    }
  }

}


void node::addEdge(int s, int f, double _w, int _k, string sI, string sF){ 
  
  edge enew = edge(s, f, _w, _k, sI, sF);

  if(s != f){
    edges.push_back(enew); 
  } else {
    if(!selfInteraction){
      edges.push_back(enew); 
      selfInteraction = true;
    }
  }

}


void node::removeEdge(int edgeKey){ 

  for(int i=0; i<edges.size(); i++){
    if( edges[i].key == edgeKey )
      edges.erase(edges.begin()+i);
  }
}

void node::updateNodeCommunity(int newCom){ c = newCom; }

inline vector<edge> node::getEdges(){ return edges; }

inline int node::getDegree(){ 

  if( degree == 0 ){
    degree = edges.size();
    return degree;
  } 
    
  return degree;

}

void node::setDegree( int newDegree ){ degree = newDegree; }

void node::print(){ cout << "we " << w << " dist " << d << " community " << c << " Degree " << getDegree() 
			 << " visited " << v << " key " << k << " (" << sID << ")" << endl; }

void node::print(fstream* fout ){ (*fout) << "we " << w  << " dist " << d << " community " << c << " Degree " 
                                          << getDegree() << " visited " << v << " key " << k << " (" << sID 
					  << ")" << endl; }

void node::printEdges(){ 
  cout << "[ "; 
  for(int i=0; i < edges.size(); i++){
    edges[i].print(); 
  }
  cout << " ]" << endl;
}

inline node node::Clone(){ return *this; }

#endif
