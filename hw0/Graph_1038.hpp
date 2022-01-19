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
class Graph {
 private:
  //vectors that contain unique ID's.
  //When we remove nodes/edges later, we can shift proxy node IDs down
	//but keep these the same.
  //Position in vec will be node position in graph, and value in vector will
	//have unique ID so we can slice into vector of proxynodes/edges
	//(will never delete element from proxynode/edge vecs)
  std::vector<int> nIDconv;//conv for conversion (from proxy id to unique id)
  std::vector<int> eIDconv;
  //Forward declaration of structs used by Node and Edge classes
  struct ProxyNodes;
  struct ProxyEdges;
  struct AdjacentNodes;//used to lessen search in has_edge(a,b))
  //vectors of proxy nodes and proxy edges
  std::vector<ProxyNodes> pnvec;
  std::vector<ProxyEdges> pedgevec;
  //total # of nodes and edges in graph
  //use unsigned to follow example given in class
  unsigned size_;
  unsigned totalEdges_;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

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
  Graph(){
      // total # of nodes and edges are 0 to start off with
      size_=0;
      totalEdges_=0;
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
  class Node {
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
    Node() {
      //already creates invalid node
    }

    /** Return this node's position. */
    const Point& position() const {
      //in this node's graph slice into correct proxy node & return Point att
      //arrow syntax can be used to access elements of a struct:
      //https://www.tutorialspoint.com/What-is-arrow-operator-in-Cplusplus
      return mygraph_->pnvec[nodeIndex_].P;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      //Same idea as position(), but return index attribute
      return mygraph_->pnvec[nodeIndex_].index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      //if index of this node and node passed in same AND graphs are same
      if(n.nodeIndex_==nodeIndex_ && mygraph_==n.mygraph_){
	return true;
	}else{
          (void) n;          // Quiet compiler warning
          return false;
	  }
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
      //if this node has index less than index of node passed in
      //will likely decide to change this as we learn more ab how < is used
      if(nodeIndex_<n.nodeIndex_){
        return true;
	}else{
           (void) n;           // Quiet compiler warning
            return false;
	  }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    //Graph needs a way to construct valid Node objects
    Graph* mygraph_;
    size_type nodeIndex_; //node id for this node
    //since we return Node() instances sometimes (in the add_node method for ex)
    //was throwing errors until I used const cast
    //https://en.cppreference.com/w/cpp/language/const_cast
    Node(const graph_type* graph, size_type index)
	:mygraph_(const_cast<graph_type*>(graph)),nodeIndex_(index) {}
  }; //END OF NODE CLASS

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
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
  Node add_node(const Point& position) {
    //initialize vector of adjacent node instances as att of proxynode struct
    std::vector<AdjacentNodes> theseAdj {};
    //add proxynode instance to graph attribute that is vector of proxynodes
    pnvec.push_back(ProxyNodes(position,size_,theseAdj));
    //add id to vec that will serve as mapping from proxyid to unique id
    //won't delete from proxynode vec so using size -1 will always be unique
    nIDconv.push_back(pnvec.size()-1);
    (void) position;      // Quiet compiler warning
    size_++;
    //return the newly created node using the private constructor.
    return Node(this, pnvec.size()-1);//use "this" to keep same graph
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //return node with the same id as the one being checked
    Node checkNode = node(n.index());
    //see if these instances are the same
    if(n==checkNode){
      return true;
      }else{
        (void) n;            // Quiet compiler warning
        return false;
        }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if(i<size()){//if index does not exceed size of graph. (calls size method)
      //use Node constructor & vector of ID's that map proxyID to a uniqueID
      return Node(this,nIDconv[i]);
    }else{
      //if the index exceeds the graph size, return an invalid node
      (void) i;             // Quiet compiler warning
      return Node();//invalid node
      }
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      //already invalid edge technically
    }

    /** Return a node of this Edge */
    Node node1() const {
      //use node1 attribute of this edge
      return Node(mygraph_,nodeClassIdx1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //use node2 attribute of this edge
      return Node(mygraph_,nodeClassIdx2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //can be equal in either "direction"
      if(e.nodeClassIdx1_ == nodeClassIdx1_ && e.nodeClassIdx2_ ==\
		 nodeClassIdx2_){
         return true;
      }else if(e.nodeClassIdx1_ == nodeClassIdx2_ && e.nodeClassIdx2_ ==\
		 nodeClassIdx1_){
         return true;
      }else{
         (void) e;           // Quiet compiler warning
         return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //just like nodes, compare if this edge id less than one we pass in
      if(edgeClassIndex_<e.edgeClassIndex_){
         return true;
      }else{
        (void) e;           // Quiet compiler warning
      return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    //Graph needs a way to construct valid Edge objects
    graph_type* mygraph_;
    //this edge has its own index, and id's for both nodes
    size_type edgeClassIndex_;
    size_type nodeClassIdx1_;
    size_type nodeClassIdx2_;
    //just like in Node constructor, was receiving error until const cast
    Edge(const graph_type* graph,size_type idx,size_type n1,size_type n2)
	:mygraph_(const_cast<graph_type*>(graph)), edgeClassIndex_(idx),\
	 nodeClassIdx1_(n1),nodeClassIdx2_(n2){}
  };//END OF EDGE CLASS

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return totalEdges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
      (void) i;             // Quiet compiler warning
      //use vector w/unique id's to define this edge's id
      return Edge(this,eIDconv[i],pedgevec[i].nodeIdx1,pedgevec[i].nodeIdx2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //obtain the index of one node
    size_type  nodeAidx=a.nodeIndex_;
     //each proxy node has vector attribute describing connected nodes
     //only search through this vector for the other node index
     for(unsigned int i = 0; i<pnvec[nodeAidx].ConnectedNodes.size(); i++){
    	if(pnvec[nodeAidx].ConnectedNodes[i].connectedNode==b.nodeIndex_)
		return true;
        }
   return false;
   (void) a; (void) b;   // Quiet compiler warning
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
      //check if edge is already in the graph
      if(has_edge(a,b)==false){
        //add the edge to graph's vector of proxy edges
        pedgevec.push_back(ProxyEdges(pedgevec.size()-1,a.nodeIndex_,\
	b.nodeIndex_));
        //add this id to the vector of id's that will map to proxy node id's
        eIDconv.push_back(pedgevec.size()-1);
        //for proxy node a, add b as a connected node
        pnvec[a.nodeIndex_].ConnectedNodes.push_back(AdjacentNodes(\
	pedgevec.size()-1,b.nodeIndex_));
        //for proxy node b, add a as a connected node
        pnvec[b.nodeIndex_].ConnectedNodes.push_back(AdjacentNodes(\
	pedgevec.size()-1,a.nodeIndex_));
        totalEdges_++;
      }
    //if edge exists in graph already, dont inc total num of edges attribute
    (void) a, (void) b;   // Quiet compiler warning
    //either way, return an Edge instance
    return Edge(this,pedgevec.size()-1,a.nodeIndex_,b.nodeIndex_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    //reset graph attributes to 0
    size_=0;
    totalEdges_=0;
    //credit to https://www.geeksforgeeks.org/vector-erase-and-clear-in-cpp/
    pedgevec.clear();
    pnvec.clear();
    nIDconv.clear();
    eIDconv.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:
   //used within an attribute of proxy nodes
   struct AdjacentNodes{
       size_type edgeID;//edge connecting the node
       size_type connectedNode;//adjacent node
       //make constructor for when we add an instance to vec in add_edge
       AdjacentNodes(size_type edgeid, size_type nextNode)
         :edgeID(edgeid),connectedNode(nextNode){}
   };

  //vector of ProxyNodes instances is graph class attribute
  struct ProxyNodes{
      Point P; // struct defined in other .hpp file
      size_type index; //proxy id. will be shifted if nodes removed
      std::vector<AdjacentNodes> ConnectedNodes;//each node has connected nodes
      //make constructor for when we add instance to vector
      //Point can be const as no need to change it
      ProxyNodes(const Point& MyP, size_type myIndex, \
	std::vector<AdjacentNodes> adjNodes)
	:P(MyP),index(myIndex), ConnectedNodes(adjNodes){}
  };

  //use same idea as for nodes since it worked
  struct ProxyEdges{
      //proxy edge id will be shifted if edges removed
      size_type edgeIndex;
      //edges have 2 nodes associated with them
      size_type nodeIdx1;
      size_type nodeIdx2;
      //make constructor for when we add to vec of proxy edges
      ProxyEdges(size_type eIndex,size_type nIdx1,size_type nIdx2)
	:edgeIndex(eIndex),nodeIdx1(nIdx1),nodeIdx2(nIdx2){}
    };

};//end of graph class

#endif // CME212_GRAPH_HPP

