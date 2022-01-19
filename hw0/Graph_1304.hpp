#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <assert.h>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

using namespace std;


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
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

  //do a map

  //(HW0 comment): constructing a struct storing info of nodes
  //include its position, index, 
  // and a vector of indexes of other nodes it is connected to (edges) 
  struct node_details{
    Point the_point;
    size_type the_index;
    std::vector<size_type> pair_indexes; //indexes of other nodes that
                                        //are connected to this node
    //default constructor with proxy design
    node_details(): the_point(), the_index(-1), pair_indexes(){}

    //constructor:
    node_details(const Point &position, size_type index):
    the_point(position), the_index(index), pair_indexes(){}

  };

  //(HW0 comment): construct a struct storing info of edges
  //include its index, the index of the two nodes it connected,
  //so we can quickly connect edge to the two nodes. 
  struct edge_details{
    size_type edge_index;
    size_type a_index;
    size_type b_index;
    //default constructor
    edge_details(): edge_index(), a_index(), b_index(){}
    //constructor with params.
    edge_details(size_type index, size_type index1, size_type index2): 
      edge_index(index), a_index(index1), b_index(index2){}
  };


  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    num_of_egdes = 0;
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
      // HW0: YOUR CODE HERE
    }

    Node(const Graph *new_graph, size_type idx):
    the_graph(const_cast<Graph*> (new_graph)), curr_index(idx) {}


    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE 
      assert(curr_index < the_graph->num_nodes()); //curr_index within range; 
      return the_graph->nodes_map[curr_index].the_point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(curr_index < the_graph->num_nodes()); //curr_index within range; 
      return curr_index;
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
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      if(this->the_graph == n.the_graph 
         && this->curr_index == n.curr_index)
      {
        return true;
      }
      else
      {
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
      // HW0: YOUR CODE HERE
      (void) n;           // Quiet compiler warning
      if(this->curr_index < n.curr_index)
      {
        return true;
      }
      else if (this->curr_index > n.curr_index)
      {
        return false;
      }
      else
      {
        if(this->the_graph < n.the_graph)
        {
          return true;
        }

        return false;

      }

    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* the_graph;
    size_type curr_index;

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_map.size();
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
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning
    nodes_map[num_nodes()] = node_details(position, num_nodes());
    return Node(this, num_nodes() - 1);  
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    if(n.the_graph == this)
    {
      return true;
    }
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
    (void) i;             // Quiet compiler warning
    return Node(this, i);     
  }

  //
  // EDGES
  //

  /** @class Graph::E dge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    Edge(const Graph* new_graph, 
         size_type new_index1, size_type new_index2):
    the_graph(const_cast<Graph*>(new_graph)), 
    the_index_1(new_index1), the_index_2(new_index2){};

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(the_graph, the_index_1);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(the_graph, the_index_2);  
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      //check graph:
      if(the_graph == e.the_graph)
      {
        if(the_index_1 == e.the_index_1 &&
           the_index_2 == e.the_index_2)
        {
          return true;
        }
        else if (the_index_1 == e.the_index_2 &&
                 the_index_2 == e.the_index_1)
        { 
          return true;
        }
        
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      if(the_index_1 + the_index_2 < e.the_index_1 + e.the_index_2)
      {
        return true;
      }
      else
      {
        return false;
      } 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* the_graph;
    size_type the_index_1;
    size_type the_index_2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_of_egdes;
  }

  //  Return the edge with index @a i.
  // @pre 0 <= @a i < num_edges()
  //   Complexity: No more than O(num_nodes() + num_edges()), hopefully less
  
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    size_type a = edges_map.at(i).a_index;
    size_type b = edges_map.at(i).b_index;
    
    return Edge(this, a, b);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning
    if(a.the_graph == b.the_graph
       && a.curr_index < nodes_map.size()
       && b.curr_index < nodes_map.size())
    {
      if(nodes_map.at(a.curr_index).pair_indexes.size() <
         nodes_map.at(b.curr_index).pair_indexes.size())
      {  
        for(size_type i = 0; i<nodes_map.at(a.curr_index).pair_indexes.size(); i++)
        {
          if(b.curr_index == nodes_map.at(a.curr_index).pair_indexes[i])
          {
            return true;
          }
        }
      }
      else
      {
        for(size_type i = 0; i<nodes_map.at(b.curr_index).pair_indexes.size(); i++)
        {
          if(a.curr_index == nodes_map.at(b.curr_index).pair_indexes[i])
          {
            return true;
          }
        }
      } 
    }
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
    (void) a, (void) b;   // Quiet compiler warning
    if(!has_edge(a,b))
    {
      
      nodes_map[a.curr_index].pair_indexes.push_back(b.curr_index);
      if(a.curr_index!= b.curr_index)
      {
        nodes_map[b.curr_index].pair_indexes.push_back(a.curr_index);
        edges_map[num_of_egdes] 
            = edge_details(num_of_egdes, a.curr_index, b.curr_index);
      }
      num_of_egdes ++;
      
    }
    
    return Edge(this, a.curr_index, b.curr_index);       
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_map.clear(); //clear maps
    edges_map.clear();
    num_of_egdes = 0;
    cout << "In Graph::clear() method: num of nodes now " 
         << num_nodes() << endl;
    cout << "        num of edges now" << num_edges() << endl;
    cout << "        num of nodes map size " << nodes_map.size() << endl;
    cout << "        num of edges map size " << edges_map.size() << endl;
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

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::map <size_type, node_details> nodes_map;
  std::map <size_type, edge_details> edges_map;
  size_type num_of_egdes;

};

#endif // CME212_GRAPH_HPP
