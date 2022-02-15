#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V, typename E>
class Graph {
 private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Predeclare internal structures
  struct internal_node;
  struct internal_edge;
 
  std::vector<internal_node> node_elements; // vector of a internal nodes
  std::vector<internal_edge> edge_elements; // vector of a internal nodes

 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  typedef V node_value_type;
  typedef E edge_value_type;

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    last_nid = 0;
    last_eid = 0;
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */

    // HW0: YOUR CODE HERE
    Node(): parent(nullptr), nid{} {};
    
    ~Node() = default; //destructor

    /** Return this node's position. */
    Point& position() const {
      // HW0: YOUR CODE HERE
      return (parent->node_elements[index()]).point;
    }

    // return nodes position, modifiable
    Point& position() {return (parent->node_elements[index()]).point;};

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      if (parent==nullptr) return size_type(-1);
      auto idx = (*parent).findnodeidx.find(nid);
      if (idx != (*parent).findnodeidx.end()) return idx->second;
      return size_type(-1);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){return parent->node_elements[index()].val;};
    
    const node_value_type& value() const{
      return parent->edge_elements[index()].val;
    };

    size_type degree() const{
      auto tonodes = parent->adjacentnodes.find(nid);
      if (tonodes == parent->adjacentnodes.end()){
        return 0;
      }else{
        return (tonodes->second).size();
      };
    };
    incident_iterator edge_begin() const {return incident_iterator(parent,nid);};
    incident_iterator edge_end() const {return incident_iterator();};

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (index() == n.index() && parent == n.parent) return true;
      return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (parent != n.parent) return (parent < n.parent);
      if (index() < n.index()) return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* parent;
    size_type nid; //unique node id. does not change.

    // this point is always const, so to call the node from within the graph, 
    // we need to define a constructor with const graph*. We use const_cast to add node or remove. 
    Node(const Graph* parent, size_type nid): parent(const_cast<Graph*>(parent)), nid(nid){};

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_elements.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type n = num_nodes();
    //internal_node* new_internalnode = new internal_node(position, n, last_nid);
    internal_node new_internalnode(position, n, last_nid, val);
    Node new_node(this,last_nid);
    findnodeidx[last_nid] = n;
    node_elements.push_back(new_internalnode);
    ++ last_nid;
    return new_node;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    Node tempnode(this,node_elements[n.index()].nid);
    if (n == tempnode) return true; 
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type n = num_nodes();
    if (i<n) return Node(this, node_elements[i].nid);
    return Node();        // Invalid node
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    // HW0: YOUR CODE HERE
    Edge(): parent(nullptr), node1_nid{}, node2_nid{}, eid{} {}

    Edge(const Graph* parent, size_type node1_nid, size_type node2_nid, size_type eid):
      parent(const_cast<Graph*>(parent)), node1_nid(node1_nid), node2_nid(node2_nid), eid(eid) {};


    ~Edge() = default; //destructor


    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(parent, node1_nid);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(parent, node2_nid);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (parent != e.parent) return false;
      if ((node1() == e.node1()) && (node2() == e.node2())) return true;
      if ((node1() == e.node2()) && (node2() == e.node1())) return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (parent != e.parent) return parent < e.parent;
      if (eid < e.eid) return true;
      return false;
    }

    // return edge index
    size_type index() const {
      if (parent==nullptr) return size_type(-1);
      auto idx = (*parent).findedgeidx.find(eid);
      if (idx != (*parent).findedgeidx.end()) return idx->second;
      return size_type(-1);
    }

    edge_value_type& value() {return parent->edge_elements[index()].val;};
    const edge_value_type& value() const {return parent->edge_elements[index()].val;};

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* parent;
    size_type node1_nid;
    size_type node2_nid;
    size_type eid;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return findedgeidx.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < num_edges()) {
      internal_edge e = edge_elements[i];
      return Edge(this, e.node1_nid, e.node2_nid, e.eid);
    }
    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    auto a1 = findedge_fromnode.find(std::make_pair(a.nid, b.nid));
    if (a1 != findedge_fromnode.end()) return true;

    auto b1 = findedge_fromnode.find(std::make_pair(b.nid, a.nid));
    if (b1 != findedge_fromnode.end()) return true;

    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    // if the graph has the node already, return the edge
    auto a1 = findedge_fromnode.find(std::make_pair(a.nid, b.nid));
    if (a1 != findedge_fromnode.end()){
      return Edge(this, a.nid, b.nid, a1->second);
    };

    auto b1 = findedge_fromnode.find(std::make_pair(b.nid, a.nid));
    if (b1 != findedge_fromnode.end()){
      return Edge(this, a.nid, b.nid, b1->second);// Invalid Edge
    };

    // otherwise, the graph does not have this edge.
    size_type m = num_edges();

    //// update edge
    internal_edge ie(a,b,m,last_eid);
    edge_elements.push_back(ie);

    //// to find edge from pair of nodes
    findedgeidx[last_eid] = m;
    findedge_fromnode[std::make_pair(a.nid,b.nid)] = last_eid;

    //// to find incident edges from the nodes
    auto vec_a = adjacentnodes.find(a.nid);
    if (vec_a == adjacentnodes.end()){
      // this node doesnt have any connected nodes
      //adjacentnodes[a.nid] = new std::vector<size_type>; //allocate on heap
      std::vector<size_type> newlist; 
      adjacentnodes[a.nid] = newlist;
      adjacentnodes[a.nid].push_back(b.nid);
    }else{
      (vec_a->second).push_back(b.nid);
    };

    auto vec_b = adjacentnodes.find(b.nid);
    if (vec_b == adjacentnodes.end()){
      // this node doesnt have any connected nodes
      //adjacentnodes[b.nid] = new std::vector<size_type> //allocate on heap
      std::vector<size_type> newlist; 
      adjacentnodes[b.nid] = newlist;
      adjacentnodes[b.nid].push_back(a.nid);
    }else{
      (vec_b->second).push_back(a.nid);
    };
    
    // return new edge 
    return Edge(this, a.nid, b.nid, last_eid++); // increment after returning
  };


  Edge FindEdge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    // if the graph has the node already, return the edge
    auto a1 = findedge_fromnode.find(std::make_pair(a.nid, b.nid));
    if (a1 != findedge_fromnode.end()){
      return Edge(this, a.nid, b.nid, a1->second);
    };

    auto b1 = findedge_fromnode.find(std::make_pair(b.nid, a.nid));
    if (b1 != findedge_fromnode.end()){
      return Edge(this, a.nid, b.nid, b1->second);// Invalid Edge
    };;

    return Edge();
  };


  /** remove_node
   * @brief removes a node and the associated data with it in the graph
   * Complexity: O(num_nodes())
   * @post g.node(i).index() == i for all i with 0 <= i < g.num_nodes()
   * @post g.node(n.index()) == n ; for all nodes not deleted
   * @post all elements from node_idx, graph_nodes and node2node in the graph that are associated with the node is destroyed
   * @post all edges from the node is removed. 
   * 
  */
  size_type remove_node(const Node&n){
    size_type idx = n.index(); 

    // if the edge does not exist, return 0
    if (n.parent==nullptr) return 0;
    if (!has_node(n)) return 0; 

    // remove associated edges
    auto connectednodes = adjacentnodes.find(n.nid);
    if (connectednodes != adjacentnodes.end()){
      std::vector<size_type> cn_vec = connectednodes->second; // copy here
      for(unsigned iijj = cn_vec.size()-1; iijj >=0 && iijj < cn_vec.size(); --iijj){
        Node n2(this,cn_vec[iijj]);
        remove_edge(n,n2);
      };
      adjacentnodes.erase(connectednodes);
    };

    // delete internal nodes
    node_elements.erase(node_elements.begin()+idx);

    // delete find nodeidx and reinindex
    findnodeidx.erase(n.nid);
    unsigned nodeidx_ne = 0; 
    for(auto it = node_elements.begin(); it != node_elements.end(); ++it){
      findnodeidx[(*it).nid] = nodeidx_ne++; //makesure update is in order
    };

    return 1;
  };

  node_iterator remove_node(node_iterator n_it){
    node_iterator nitcopy(n_it.parent, n_it.thisnode);
    nitcopy++;
    size_type sz = remove_node(*n_it);
    return nitcopy;
  };

    /** remove_edge
   * @brief removes an edge and the associated data with it in the graph
   * Complexity: O(num_nodes() + num_edges())
   * @post num_edges() returns number of unique undirected edges
   * 
  */
  size_type remove_edge(const Edge & e){
    //useful data
    Node n1 = e.node1();
    Node n2 = e.node2();
    size_type idx = e.index(); 

    // if the edge does not exist, return 0
    if (e.parent==nullptr) return 0;
    if (!has_edge(n1,n2)) return 0; 

    // delete associated internal_edge
    edge_elements.erase(edge_elements.begin()+idx);
    findedgeidx.erase(e.eid);

    // iterate through edge_elements and reindex findedgeidx map O(num_edges)
    unsigned edgeidx_ee = 0; 
    for(auto it = edge_elements.begin(); it != edge_elements.end(); ++it){
      findedgeidx[(*it).eid] = edgeidx_ee++; //makesure update is in order
    };

    // iterate through findedge_fromnode O(num_nodes)
    auto a1 = findedge_fromnode.find(std::make_pair(n1.nid, n2.nid));
    if (a1 != findedge_fromnode.end()) findedge_fromnode.erase(a1);

    auto b1 = findedge_fromnode.find(std::make_pair(n2.nid, n1.nid));
    if (b1 != findedge_fromnode.end()) findedge_fromnode.erase(b1);

    // remove adjacent nodes 
    std::vector<size_type>& nodesvec1 = adjacentnodes[n1.nid];
    for (auto an_iter = nodesvec1.begin();an_iter!=nodesvec1.end();an_iter++){
      if (*an_iter == n2.nid){
        nodesvec1.erase(an_iter); //erase node1 from the vector of nodes connected to node2
        break;
      };
    };

    std::vector<size_type>& nodesvec2 = adjacentnodes[n2.nid];
    for (auto an_iter = nodesvec2.begin();an_iter!=nodesvec2.end();an_iter++){
      if (*an_iter == n1.nid){
        nodesvec2.erase(an_iter); //erase node1 from the vector of nodes connected to node2
        break;
      };
    };

    return 1;
  };

  size_type remove_edge(const Node& n1, const Node& n2){
    return remove_edge(FindEdge(n1,n2));
  };

  edge_iterator remove_edge(edge_iterator e_it){
    edge_iterator eitcopy(e_it.parent, e_it.currEdge);
    eitcopy++;
    size_type sz = remove_edge(*e_it);
    return eitcopy;
  };

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_elements.clear(); node_elements = {};
    edge_elements.clear(); edge_elements = {};
    findnodeidx ={};
    findedgeidx = {};
    findedge_fromnode = {}; 
    last_eid = 0;
    last_nid = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator(): parent(nullptr), thisnode(Node()) {}; 

      
    // constructor
    NodeIterator(const Graph* graph, value_type somenode): 
      parent(const_cast<Graph*>(graph)), thisnode(somenode) {};


    // HW1 #2: YOUR CODE HERE
    Node operator*() const {return thisnode;}; 

    NodeIterator& operator++() {
      // return the next node in index
      size_type idx = (thisnode.index()) + 1;
      if (idx < (*parent).num_nodes()){
        Node nextnode(parent, parent->node_elements[idx].nid);
        thisnode = nextnode;
      }else{
        thisnode = Node();
        parent = nullptr;
      };
      return *this;

    };

    bool operator==(const NodeIterator& otherNI) const {
      return thisnode == otherNI.thisnode;
    };

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* parent;        
    value_type thisnode;  
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    if (size()>0){
      Node firstnode(this, node_elements[0].nid);
      return NodeIterator(this,firstnode);
    }else{
      return NodeIterator();
    };
  };

  node_iterator node_end() const {return NodeIterator();};

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator(): parent(nullptr), source_nid(0), idx(0), curr_edge(Edge()) {};

    IncidentIterator(const Graph* initpar, size_type nid): 
      parent(nullptr), source_nid(0), idx(0), curr_edge(Edge()){
        Node fromnode(initpar, source_nid);
        if (fromnode.degree() > 0){
          parent = const_cast<Graph*>(initpar);
          source_nid = nid;
          curr_edge = findEdge(idx);
        };
      };

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {return curr_edge;};

    Edge findEdge(size_type idx){
      auto next_nids = parent->adjacentnodes.find(source_nid);
      if (next_nids == parent->adjacentnodes.end()){
        return Edge();
      };
      size_type nid = (next_nids->second)[idx]; 
      size_type eid; 
      Edge e;
      auto e_a = parent->findedge_fromnode.find(std::make_pair(source_nid, nid));
      if (e_a != parent->findedge_fromnode.end()){
        eid = e_a->second;        
        e = Edge(parent, source_nid, nid, eid);
      }else{
        auto e_b = parent->findedge_fromnode.find(std::make_pair(nid, source_nid)); // a little inefficient to search both..
        if (e_b != parent->findedge_fromnode.end()){
          eid = e_b->second;
          e = Edge(parent, source_nid, nid, eid);
        }else{
          e = Edge(); // should not return this if there is an edge.
        };
      };
      return e;
    };
    
    IncidentIterator& operator++() {
      if (curr_edge.parent == nullptr) {return *this;};
      
      idx++;
      Node fromnode(parent, source_nid);
      size_type deg = fromnode.degree();
      if ((idx < deg) && (idx > 0)){
        auto tonodes = parent->adjacentnodes.find(source_nid);
        if (tonodes == parent->adjacentnodes.end()){
          return *this;
        };
        curr_edge = findEdge(idx);        
      }else{
        parent = nullptr;
        source_nid = 0;
        idx = 0;
        curr_edge = Edge();
      };
      
      return *this;
    };
    
    bool operator==(const IncidentIterator& someIncIter) const {
      return ((parent == someIncIter.parent) &&
              (source_nid == someIncIter.source_nid) && 
              (idx == someIncIter.idx) &&
              (curr_edge == someIncIter.curr_edge));
    };

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* parent;
    size_type source_nid; //nid of the source node.
    size_type idx;
    value_type curr_edge;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(): parent(nullptr), currEdge(Edge()) {};

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {return currEdge;};
    EdgeIterator& operator++(){
      if (parent == nullptr) {return *this;};

      size_type idx = currEdge.index() + 1;
      size_type num_edges = parent->num_edges();
      if (idx < num_edges){
        internal_edge ie = (parent->edge_elements[idx]);
        currEdge = Edge(parent,ie.node1_nid, ie.node2_nid,ie.eid);
      }else{
        // *this = EdgeIterator();
        parent = nullptr;
        currEdge = Edge();
        idx = 0;
      };
      return *this;

    };

    bool operator==(const EdgeIterator& someEdgeIter) const {
      return ((parent == someEdgeIter.parent)&&
              (currEdge == someEdgeIter.currEdge));
    };

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* parent;
    value_type currEdge;

    EdgeIterator(const Graph* parent, value_type init_edge): 
      parent(const_cast<Graph*>(parent)), currEdge(init_edge) {};
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    if (edge_elements.size() < 1){
      return edge_iterator();
    }
    internal_edge e =  edge_elements[0];
    Edge firstedge(this, e.node1_nid,  e.node2_nid, e.eid);
    return edge_iterator(this, firstedge);
  };

  edge_iterator edge_end() const {return edge_iterator();};

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node{
    Point point;
    size_type node_idx;
    size_type nid;
    node_value_type val;

    // constructor without val.
    internal_node(const Point& point, size_type node_idx, size_type nid) 
      : point(point), node_idx(node_idx), nid(nid), val() {};

    internal_node(const Point& point, size_type node_idx, size_type nid, node_value_type val) 
    : point(point), node_idx(node_idx), nid(nid), val(val) {};

  };

  struct internal_edge{
    size_type node1_nid;
    size_type node2_nid;
    size_type edge_idx;
    size_type eid;
    edge_value_type val; 

    internal_edge(const Node& node1, const Node& node2, size_type edge_idx, size_type eid) : 
      node1_nid(node1.nid), node2_nid(node2.nid), edge_idx(edge_idx), eid(eid), val() {};

    internal_edge(const Node& node1, const Node& node2, size_type edge_idx, size_type eid, edge_value_type val) : 
      node1_nid(node1.nid), node2_nid(node2.nid), edge_idx(edge_idx), eid(eid), val(val) {};
  };

  std::map<size_type, size_type> findnodeidx; // unique nid to index in graph
  std::map<size_type, size_type> findedgeidx; // unique eid to index in graph
  std::map<std::pair<size_type,size_type>, size_type> findedge_fromnode;

  std::map<size_type, std::vector<size_type>> adjacentnodes; // unique nid to vector of nid its connected by an edge

  size_type last_nid;
  size_type last_eid;

};

#endif // CME212_GRAPH_HPP
